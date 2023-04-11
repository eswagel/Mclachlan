include("HermiteSymbolic.jl");
include("../JuliaMathematica/SymbolicsMathLink.jl");
using DifferentialEquations
import DASKR: daskr
import Roots: find_zero
#using IntervalArithmetic, IntervalRootFinding
#import Sundials: IDA
#using Plots

genPsi(power::Int, dn::Vector{Num},skip::Int=1)=begin
    ret = Vector{Terms}(undef,power+1)
    for i=0:power
        ret[i+1] = Terms(Hermite(1,-dn[i+1],0,i*skip))#hermite_poly(skip*i,X)*Hermite(1, -dn[i+1], 0, 0)
    end
    ret
end

function calcKsol(MclachlanX::Terms, psi::Vector{Hermite},power::Integer,K::Vector{Num},t::Num,Lop::Operator)

    #Mclachlan 2.4
    Kequations::Vector{Equation}=[GaussianIntegral(MclachlanX*psi[i])~0 for i=1:power+1]

    #Solve them with initial condition of even weights
    #This is the step that freezes but that Mathematica can do
    initconds = [1;zeros(Int,length(K)-1)]#fill((1/(power+1))::Float64,length(K))#
    DSolveMathematica(Kequations, K, t; ics=(substitute.(K,Ref(t=>0))::Vector{Num}).~initconds)
end;

subsK(power::Int, Ksol::Vector{Num},dn::Vector{Num},skip::Int=1)::Terms=Terms([Ksol[i]*Hermite(1,-dn[i],0,skip*(i-1)) for i=1:power+1])

function performMclachlan(power::Int, Lop::Operator, params::Vector{Num}, skip::Int=1)
    @variables t Karr(t)[1:power+1] dtarr(t)[1:power+1] dKarr[1:power+1] ddtarr(t)[1:power+1]
    K::Vector{Num}, dt::Vector{Num}, dK::Vector{Num}, ddt::Vector{Num} = scalarize.((Karr, dtarr, dKarr, ddtarr)::NTuple{4,Symbolics.Arr{Num, 1}})::NTuple{4, Vector{Num}}

    #psi::Vector{Hermite}=@time "genPsi" genPsi(power,dn,useEvens)
    psivector::Vector{Terms}=genPsi(power,dt,skip)
    psi::Terms=sum([psivector[i]*K[i] for i=1:power+1])

    #substitution::Dict{Num,Num} = Dict([Differential(t).(K);Differential(t).(dt)].=>[dK;ddt])
    #reverse_substitution::Dict{Num,Num} = Dict([dK;ddt].=>[Differential(t).(K);Differential(t).(dt)])

    Theta::Terms=derivative(psi, t) #in dt form
    #Thetan=substitute(Theta,dtTodn) #in n form
    MclachlanX::Terms = Theta-Lop(psi)
    #MclachlanXn=substitute(MclachlanX,dtTodn)
    
    #Mclachlan 2.4
    Kequations::Vector{Num}=[GaussianIntegral(MclachlanX*psivector[i]) for i=1:power+1]

    #Now substitute in our solutions for the weights
    #psisubs::Terms=subsK(power,Ksol,dn,useEvens)

    #println("Substituted")
    #Setting up Mclachlan 2.5
    #substitution = Dict(K.=>Ksol)
    #Mclachlan=@time "genMclachlan" genMclachlan(dn,(@time "substituteTheta" substitute(Theta, substitution)),(@time "substituteMclachlanX" substitute(MclachlanX,substitution)))

    #Mclachlan 2.5
    mclachlaneqs = @. GaussianIntegral(derivative($Ref(Theta),dt) * $Ref(MclachlanX),$Val(true))
    return (expr=[Kequations;mclachlaneqs], dt=dt, K=K, dK=dK, ddt=ddt)
end;

function calcu0(H::Operator)::Vector{Float64}
    @variables dn
    mdn = W"dn"
    psi = Hermite(1, -dn, 0, 0)

    energy = GaussianIntegral(psi*(H*psi))

    result::Vector{Pair{Num,Num}} = wcall("FindRoot", derivative(energy,dn), [dn,1.0,0.0,"Infinity"])
    Float64[result[i][2].val::Float64 for i=1:length(result)]
end;

function solveMclachlanForSteadyState(power::Int, K::Vector{Num},dt::Vector{Num},Hsubs::Operator,resultfunc1::Function,param_vals::Vector{Float64}, variationald0::Float64)
    findrootinit::Vector{Num}=resultfunc1(fill(0.0,2(power+1)),[K;dt],param_vals,0.0)

    initcondsguess=[[[K[1]::Num,1.0]];[[K[i]::Num,1.0,0.0,W"Infinity"] for i=2:power+1];[[dt[i]::Num,(i>1 ? 1.0 : 1.0),0.0,W"Infinity"] for i=1:power+1]]
    solpairsarray=wcall("FindRoot",findrootinit,initcondsguess,MaxIterations=1000)::Vector{Pair{Num,Num}}
    Float64[real(solpairsarray[i][2]).val for i=1:length(solpairsarray)]
end

function zeroDerivativesSolve(tfinal::Float64,power::Int, skip::Int, K::Vector{Num},dt::Vector{Num},H::Operator,resultfunc::Tuple{Function,Function},param_vals::Vector{Float64}, variationald0::Float64)
    param_subs = Dict(params::Vector{Num}.=>param_vals)
    Hsubs = substitute(H,param_subs)
    println(param_vals)

    steadystate=solveMclachlanForSteadyState(power, K, dt, Hsubs, resultfunc[1], param_vals, variationald0)
    ksol=steadystate[1:power+1]
    dsol=steadystate[power+2:end]

    return (ksol, dsol)
end

function initcondsSolveMclachlan(power::Int, K::Vector{Num},dt::Vector{Num},Hsubs::Operator,resultfunc1::Function,param_vals::Vector{Float64}, variationald0::Float64)
    findrootinit::Vector{Num}=resultfunc1([0.1;fill(0.1,2(power+1)-1)],[K;dt],param_vals,0.0)

    initcondsguess=[[[K[1]::Num,1.0,0.0,"Infinity"]];[[K[i]::Num,1.0,0.0,"Infinity"] for i=2:power+1];[[dt[i]::Num,(i>1 ? variationald0 : variationald0),0.0,"Infinity"] for i=1:power+1]]
    solpairsarray=wcall("FindRoot",findrootinit,initcondsguess,MaxIterations=100)::Vector{Pair{Num,Num}}
    Float64[solpairsarray[i][2].val for i=1:length(solpairsarray)]
end

function flowingTimeSolve(tfinal::Float64,power::Int, skip::Int, K::Vector{Num},dt::Vector{Num},H::Operator,resultfunc::Tuple{Function,Function},param_vals::Vector{Float64}, variationald0::Float64)
    param_subs = Dict(params::Vector{Num}.=>param_vals)
    Hsubs = substitute(H,param_subs)

    
    initconds = initcondsSolveMclachlan(power, K, dt, Hsubs, resultfunc[1], param_vals, variationald0)
    println(param_vals)

    daef=DAEFunction(resultfunc[2],syms=[:K,:d],paramsyms=[:a,:c,:g])
    daep=DAEProblem(daef,[0.1;fill(0.1,2(power+1)-1)],initconds,[0.0,tfinal],param_vals)
    try
        daesol=@time "Solving DAE" solve(daep,daskr())

        ksol::Vector{Float64}=abs.(daesol[end][1:power+1])
        dsol::Vector{Float64}=daesol[end][power+2:end]

        return (ksol, dsol)
    catch e
        println("Error in DAE solve")
        println(e)
        println("Returning NaNs")
        return (fill(NaN,power+1),fill(NaN,power+1))
        #return zeroDerivativesSolve(tfinal,power, skip, K, dt, H, resultfunc, param_vals, variationald0)
    end

end

function transform_result_to_function(mclachlanResult)
    ret = Vector{Tuple{Function,Function}}(undef,length(mclachlanResult))
    for i=1:length(mclachlanResult)
        result = mclachlanResult[i]
        ret[i] = eval.(build_function(substitute(result.expr,Dict(Differential(t).([result.K;result.dt]).=>[result.dK;result.ddt]))::Vector{Num},[result.dK;result.ddt],[result.K;result.dt],params,t))::Tuple{Function,Function}
    end
    ret
end
function solveEquations(tfinal, H, power, mclachlanResults, resultfuncs, param_vals, variationald0::Float64, method, skip)
    ksol,dsol = method(tfinal, power, skip, mclachlanResults[power+1].K, mclachlanResults[power+1].dt, H, resultfuncs[power+1],param_vals,variationald0)

    param_subs = Dict(params::Vector{Num}.=>param_vals)
    Hsubs = substitute(H,param_subs)
    if power==0 && ksol[1]==0
        ksol[1]=1.0
    end
    hermitesol=Terms([Hermite(ksol[i+1],-dsol[i+1],0,skip*i) for i=0:power]);
    l1norm = GaussianIntegral(hermitesol)
    l2norm = sqrt(GaussianIntegral(hermitesol*hermitesol))
    hermitesoll1::Terms=hermitesol/l1norm;
    hermitesoll2::Terms=hermitesol/l2norm;

    #Check that the normalization is correct
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