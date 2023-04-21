include("HermiteSymbolic.jl");
include("../JuliaMathematica/SymbolicsMathLink.jl");
using DifferentialEquations
import DASKR: daskr

genPsi(power::Int, dn::Vector{Num},skip::Int=1)=begin
    """ Generate a list of power+1 Hermite basis states. """
    ret = Vector{Terms}(undef,power+1)
    for i=0:power
        ret[i+1] = Terms(Hermite(1,-dn[i+1],0,i*skip))#hermite_poly(skip*i,X)*Hermite(1, -dn[i+1], 0, 0)
    end
    ret
end


function performMclachlan(power::Int, Lop::Operator, skip::Int=1)
    """
    Generate Mclachlan's equations for the given power and operator.

    @param power The power of the Hermite basis to use.
    @param Lop The operator to use.
    @param skip The prefactor in front of each Hermite basis state is x^(skip*i) for i=0:power.
    """
    @variables t Karr(t)[1:power+1] dtarr(t)[1:power+1] dKarr[1:power+1] ddtarr(t)[1:power+1]
    #Make these Vector{Num} instead of Symbolics.Arr{Num,1}
    K::Vector{Num}, dt::Vector{Num}, dK::Vector{Num}, ddt::Vector{Num} = scalarize.((Karr, dtarr, dKarr, ddtarr)::NTuple{4,Symbolics.Arr{Num, 1}})::NTuple{4, Vector{Num}}

    psivector::Vector{Terms}=genPsi(power,dt,skip)
    psi::Terms=sum([psivector[i]*K[i] for i=1:power+1])

    Theta::Terms=derivative(psi, t)
    MclachlanX::Terms = Theta-Lop(psi)
    
    #Mclachlan 2.4
    Kequations::Vector{Num}=[GaussianIntegral(MclachlanX*psivector[i]) for i=1:power+1]

    #Mclachlan 2.5
    mclachlaneqs = @. GaussianIntegral(derivative($Ref(Theta),dt) * $Ref(MclachlanX),$Val(true))
    return (expr=[Kequations;mclachlaneqs], dt=dt, K=K, dK=dK, ddt=ddt)
end;

function calcu0(H::Operator)::Vector{Float64}
    """ Calculate the Quantum Variational energy <H> for a Gaussian ansatz. """
    @variables dn
    mdn = W"dn"
    psi = Hermite(1, -dn, 0, 0)

    energy = GaussianIntegral(psi*(H*psi))

    result::Vector{Num} = wcall("SolveValues", derivative(energy,dn)~0, mdn, "Reals")
    Float64[result[i].val::Float64 for i=1:length(result)]
end;

function transform_result_to_function(mclachlanResult)
    """ Use build_function to take Mclachlan's equations from performMclachlan and turn them into a function. """
    ret = Vector{Tuple{Function,Function}}(undef,length(mclachlanResult))
    for i=1:length(mclachlanResult)
        result = mclachlanResult[i]
        ret[i] = eval.(build_function(substitute(result.expr,Dict(Differential(t).([result.K;result.dt]).=>[result.dK;result.ddt]))::Vector{Num},[result.dK;result.ddt],[result.K;result.dt],params,t))::Tuple{Function,Function}
    end
    ret
end
function solveEquations(tfinal, H::Operator, power, mclachlanResults, resultfuncs, param_vals, variationald0::Float64, method::Function, skip)
    """ Solve the Mclachlan equations using the specified method. """

    #Solve the equations using method. The method returns a tuple of (ksol,dsol), two Vector{Float64} of length power+1.
    ksol::Vector{Float64},dsol::Vector{Float64} = method(tfinal, power, skip, mclachlanResults[power+1].K, mclachlanResults[power+1].dt, H, resultfuncs[power+1],param_vals,variationald0)

    param_subs = Dict(params::Vector{Num}.=>param_vals)
    Hsubs = substitute(H,param_subs)
    hermitesol=Terms([Hermite(ksol[i+1],-dsol[i+1],0,skip*i) for i=0:power]);
    #Calculate normalization constants
    l1norm = GaussianIntegral(hermitesol)
    l2norm = sqrt(GaussianIntegral(hermitesol*hermitesol))
    hermitesoll1::Terms=hermitesol/l1norm;
    hermitesoll2::Terms=hermitesol/l2norm;

    #Check that the normalization is correct (should be 1)
    if(!(GaussianIntegral(hermitesoll1).val[1] ≈ 1.0))
        println("Normalization is not correct for order", power, ": ", param_vals)
        println(hermitesol)
    end

    energy::Float64 = GaussianIntegral(hermitesoll2*(Hsubs*hermitesoll2)).val
    variance::Float64 = GaussianIntegral(X^2 * hermitesoll1).val
    kurtosis::Float64 = GaussianIntegral(X^4 * hermitesoll1).val

    return (ksol=ksol./l1norm,dsol=-1 .*dsol,hermitesoll1=hermitesoll1,hermitesoll2=hermitesoll2,energy=energy,variance=variance,kurtosis=kurtosis,params=param_vals)
end
function solveEquationsParamsTable(tfinal, H, mclachlanResults, resultfuncs, paramstable, variationald0s, method, skip)
    """ Solve the Mclachlan equations using solveEquations for a table of parameters. """
    ret = Array{NamedTuple{(:ksol, :dsol, :hermitesoll1, :hermitesoll2, :energy, :variance, :kurtosis, :params), Tuple{Vector{Num}, Vector{Float64}, Terms, Terms, Float64, Float64, Float64, Vector{Float64}}}}(undef, size(paramstable)..., length(mclachlanResults))
    for i=1:size(ret,1)
        for j=1:size(ret,2)
            for k=1:length(mclachlanResults)
                ret[i,j,k] = solveEquations(tfinal, H, k-1, mclachlanResults, resultfuncs, paramstable[i,j], variationald0s[i,j], method, skip)
            end
        end
    end
    ret
end

const params = @variables a c g;
#a, c, g = (5,1,2);