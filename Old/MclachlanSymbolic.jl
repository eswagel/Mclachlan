using DifferentialEquations
import DASKR: daskr
import Plots: plot
using DataFrames: DataFrame
using CSV
using StaticArrays
using SnoopCompile
using ProfileView
using JET
using Profile

const cwd=splitdir(Base.source_path())[1]

include("HermiteSymbolic.jl");
include("../JuliaMathematica/SymbolicsMathLink.jl");

genPsi(power::Int, K::Vector{T}, dn::Vector{Num},even::Bool=false) where T= begin
    if even
        return Hermite[Hermite(K[i+1], -dn[i+1], 0, 2i) for i=0:power]
    else
        return Hermite[Hermite(K[i+1], -dn[i+1], 0, i) for i=0:power]
    end
end

function calcKsol(Kpsi::Vector{Hermite},power::Integer,K::Vector{Num},t::Num,Lop::Operator)
    #Put the L operator (it's a parameter of the function) into the basis of Hermites
    Lnm::Vector{Vector{Num}}=[[GaussianIntegral(Kpsi[n]*Lop(Kpsi[m])) for m in 1:power+1] for n in 1:power+1]
    #Mclachlan 2.4 left-hand side
    KDtK::Vector{Num}=[GaussianIntegral(derivative(Kpsi[i],t)*Kpsi[i]) for i in 1:power+1]

    #Mclachlan 2.4
    Kequations::Vector{Equation}=Equation[KDtK[i]~sum(Lnm[i])[1] for i in 1:size(Lnm)[1]]

    #Solve them with initial condition of even weights
    #This is the step that freezes but that Mathematica can do
    initconds = [1;zeros(Int,length(K)-1)]#fill((1/(power+1))::Float64,length(K))#
    DSolveMathematica(Kequations, K, t; ics=(substitute.(K,Ref(t=>0))::Vector{Num}).~initconds)
end;

subsK(power::Int, Ksol::Vector{Num},dn::Vector{Num}, dt::Vector{Num},even::Bool=false)::Terms=begin
    if even
        return Terms([substitute(Ksol[i], dn, dt)*Hermite(1,-dt[i],0,2(i-1)) for i=1:power+1])
    else
        return Terms([substitute(Ksol[i], dn, dt)*Hermite(1,-dt[i],0,(i-1)) for i=1:power+1])
    end
end;

function genMclachlan(t::Num, dn::Vector{Num}, dt::Vector{Num}, ddt::Vector{Num}, Lop::Operator, psisubs::Terms)
    Theta::Terms=derivative(psisubs,t) #In dt format
    difft::Vector{Num} = Differential(t).(dt)
    Thetaddt::Terms = substitute(Theta,difft,ddt) #In ddt format
    LpsiMinusTheta::Terms=Lop*psisubs-Theta #In dt format

    return @. GaussianIntegral(substitute(derivative($Ref(Thetaddt),dt)::Vector{Terms},$Ref(ddt),$Ref(difft))::Vector{Terms} * $Ref(LpsiMinusTheta),Val(true))
end;

function performMclachlan(power::Int, Lop::Operator, params::Vector{Num}, useEvens::Bool=false)
    @variables t Karr(t)[1:power+1] dtarr(t)[1:power+1] dnarr[1:power+1] ddtarr[1:power+1]
    K::Vector{Num}, dt::Vector{Num}, dn::Vector{Num}, ddt::Vector{Num} = scalarize.((Karr, dtarr, dnarr, ddtarr)::NTuple{4,Symbolics.Arr{Num, 1}})::NTuple{4, Vector{Num}}

    psi::Vector{Hermite}=genPsi(power,K,dn,useEvens)
    
    Ksol::Vector{Num}=calcKsol(psi,power,K,t,Lop)
    #Now substitute in our solutions for the weights
    psisubs::Terms=subsK(power,Ksol,dn,dt,useEvens)

    println("Substituted")
    #Setting up Mclachlan 2.5
    Mclachlan=genMclachlan(t,dn,dt,ddt,Lop,psisubs)

    return (Ks=Ksol, expr=Mclachlan, dt=dt, ddt=ddt)
end;

function calcu0(power::Int, H::Operator)::Vector{Num}
    @variables dn
    mdn = W"dn"
    psi = Hermite(1, -dn, 0, power)

    energy = expr_to_mathematica(GaussianIntegral(psi*(H*psi)))

    mathematica_to_expr(weval(W"SolveValues"(W"Equal"(W"D"(energy,mdn),0),mdn,W"Reals")))
end;

params = @variables a c g;

param_vals = [5.0, 1.0, 1.0]
#a, c, g = (5,1,2);

Lop=Dx*(a*X+c*X^3)+(g*Dx^2);
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;

tfinal=10.0;

result=@time performMclachlan(0, Lop, params, true);
#@code_warntype performMclachlan(1, Lop, params, true)

mdt = expr_to_mathematica(result.dt)

mexpr = @time W"Thread"(W"Equal"(expr_to_mathematica(result.expr),0));
initconds=weval(W"Table"(W"Equal"(W"Part"(W"ReplaceAll"(mdt,W`t->0`),W"i"),W"Part"(W"List"(2.56, 1.2), W"i")),W"List"(W"i",1,W"Length"(mdt))))

mexprWithParams=@time weval(mexpr; Dict(Symbol.(params).=>param_vals)...);

ndsol=@time weval(W"NDSolveValue"(W"Join"(mexprWithParams, initconds),mdt, W`{t,0,$tfinal}`,W`Method->{"EquationSimplification"->"Residual"}`));
weval(W"ReplaceAll"(ndsol, W`t->$tfinal`))
calcu0(0, substitute(Ldag*Lop, Dict(params.=>param_vals)))


(;Ks,solution)=result;
solution(tfinal)
df=DataFrame(solution)

@variables t dn[1:length(Ks)]
Kfunc=eval(Base.remove_linenums!(build_function(Ks,[t;dn])[1]))
Knumerics=(Kfunc ∘ (x->[x...])).(collect(zip(solution.t,[solution[i,:] for i in 1:length(Ks)]...)))

CSV.write(cwd*"/solution.csv",df)
CSV.write(cwd*"/Ks.csv",DataFrame(mapreduce(permutedims, vcat, Knumerics),:auto))

plot(solution)

@time sum(rand(Float64,1000))