include("HermiteSymbolic.jl");
include("../JuliaMathematica/SymbolicsMathLink.jl");
using DifferentialEquations
import DASKR: daskr
using IntervalArithmetic, IntervalRootFinding
#import Sundials: IDA
using Plots

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

function genMclachlan(dn::Vector{Num}, Theta::Terms, MclachlanX::Terms)
    return @. GaussianIntegral(derivative($Ref(Theta),dn) * $Ref(MclachlanX),$Val(true))
end;

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
    mclachlaneqs = genMclachlan(dt,Theta,MclachlanX)
    return (expr=[Kequations;mclachlaneqs], dt=dt, K=K, dK=dK, ddt=ddt)
end;

function calcu0(H::Operator)::Vector{Float64}
    @variables dn
    mdn = W"dn"
    psi = Hermite(1, -dn, 0, 0)

    energy = GaussianIntegral(psi*(H*psi))

    Symbolics.value.(wcall("SolveValues",derivative(energy,dn)~0,dn,"Reals"))
end;

function initcondsSolveMclachlan(power::Int, K::Vector{Num},dt::Vector{Num},Hsubs::Operator,resultfunc1::Function,param_vals::Vector{Float64})
    findrootinit::Vector{Num}=resultfunc1(fill(-0.01,2(power+1)),[K;dt],param_vals,0.0)
    
    variationald0=calcu0(Hsubs)[1].val

    initcondsguess=[[[K[1]::Num,1.0]];[[K[i]::Num,1.0] for i=2:power+1];[[dt[i]::Num,i>1 ? 3.0 : variationald0::Float64] for i=1:power+1]]
    solpairsarray=wcall("FindRoot",findrootinit,initcondsguess,MaxIterations=1000)::Array{Pair{Num,Float64}}
    Float64[solpairsarray[i][2] for i=1:length(solpairsarray)]
end

function solveMclachlanForSteadyState(power::Int, K::Vector{Num},dt::Vector{Num},Hsubs::Operator,resultfunc1::Function,param_vals::Vector{Float64})
    findrootinit::Vector{Num}=resultfunc1(fill(0.0,2(power+1)),[K;dt],param_vals,0.0)

    variationald0=calcu0(Hsubs)[1]

    initcondsguess=[[[K[1]::Num,1.0,1e-10,10.0]];[[K[i]::Num,1.0,0.0,W"Infinity"] for i=2:power+1];[[dt[i]::Num,(i>1 ? 3.0 : variationald0::Float64)] for i=1:power+1]]
    solpairsarray=wcall("FindRoot",findrootinit,initcondsguess,MaxIterations=1000)::Array{Pair{Num,Float64}}
    Float64[solpairsarray[i][2] for i=1:length(solpairsarray)]
end

function solveMclachlanEquations(tfinal::Float64,power::Int, skip::Int, K::Vector{Num},dt::Vector{Num},H::Operator,resultfunc1::Function,daef::DAEFunction,param_vals::Vector{Float64})
    param_subs = Dict(params::Vector{Num}.=>param_vals)
    Hsubs = substitute(H,param_subs)

    #=
    initconds = initcondsSolveMclachlan(power, K, dt, Hsubs, resultfunc1, param_vals)

    daep=DAEProblem(daef,fill(-0.01,2(power+1)),initconds,[0.0,tfinal],param_vals)
    daesol=solve(daep,daskr())

    ksol::Vector{Float64}=abs.(daesol[end][1:power+1])
    dsol::Vector{Float64}=daesol[end][power+2:end]
    =#

    steadystate=solveMclachlanForSteadyState(power, K, dt, Hsubs, resultfunc1, param_vals)
    ksol=steadystate[1:power+1]
    dsol=steadystate[power+2:end]

    hermitesol=Terms([Hermite(ksol[i+1],-dsol[i+1],0,skip*i) for i=0:power]);
    l1norm = GaussianIntegral(hermitesol)
    l2norm = sqrt(GaussianIntegral(hermitesol*hermitesol))
    hermitesoll1::Terms=hermitesol/l1norm;
    hermitesoll2::Terms=hermitesol/l2norm;

    #Check that the normalization is correct
    @assert GaussianIntegral(hermitesoll1).val[1] ≈ 1.0

    energy::Float64 = GaussianIntegral(hermitesoll2*(Hsubs*hermitesoll2)).val[1]
    variance::Float64 = GaussianIntegral(X^2 * hermitesoll1).val[1]
    kurtosis::Float64 = GaussianIntegral(X^4 * hermitesoll1).val[1]

    return (ksol=ksol./l1norm,dsol=-1 .*dsol,hermitesoll1=hermitesoll1,hermitesoll2=hermitesoll2,energy=energy,variance=variance,kurtosis=kurtosis)
end

params = @variables a c g;

param_vals = [5.0, 1.0, 1.0]
#a, c, g = (5,1,2);

Lop=Dx*(a*X+c*X^3)+(g*Dx^2);
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;
H = Ldag*Lop;

tfinal=10.0;

power,skip = 0,2
result=@time performMclachlan(power, Lop, params, skip);

resultfunc=eval.(build_function(substitute(result.expr,Dict(Differential(t).([result.K;result.dt]).=>[result.dK;result.ddt])),[result.dK;result.ddt],[result.K;result.dt],params,t))

daef=DAEFunction(resultfunc[2],syms=[:K,:d],paramsyms=[:a,:c,:g])

numij = 2;
paramtable = [[10.0, 2i+1.0,2j+1.0] for i=0:numij, j=0:numij]
solutions = [solveMclachlanEquations(tfinal, power, skip, result.K, result.dt, H, resultfunc[1], daef, [10.0, 2i+1.0,2j+1.0]) for i=0:numij, j=0:numij];

map(x->x.dsol,solutions)
map(calcu0,map(x->substitute(H,Dict(params.=>x)),paramtable))

