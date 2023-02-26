using DifferentialEquations
using DASKR
import Plots: plot
using DataFrames
using CSV

const cwd=splitdir(Base.source_path())[1]

include("Hermite.jl");
include("JuliaMathematica.jl");

function calcKsol(Kpsi::Vector{Hermite},power::Integer,K::Vector{SymPy.Sym},Lop::Operator)::Vector{SymPy.Sym}
    #Put the L operator (it's a parameter of the function) into the basis of Hermites
    Lnm::Vector{Vector{SymPy.Sym}}=[[GaussianIntegral(Kpsi[n]*Lop(Kpsi[m])) for m in 1:power+1] for n in 1:power+1]
    #Mclachlan 2.4 left-hand side
    KDtK::Vector{SymPy.Sym}=[GaussianIntegral(SymPy.diff(Kpsi[i],t)*Kpsi[i]) for i in 1:power+1]

    #Mclachlan 2.4
    Kequations::Vector{SymPy.Sym}=Sym[SymPy.Eq(KDtK[i],sum(Lnm[i])[1]::Sym)::Sym for i in 1:size(Lnm)[1]]
    #Solve them with initial condition of even weights
    #This is the step that freezes but that Mathematica can do
    DSolveMathematica(Kequations, K; ics=Dict(K[i](0)::Sym=>1/sqrt(power+1) for i in 1:power+1))
end;

function calcu0(power::Integer, dn::Vector{SymPy.Sym}, Lop::Operator, Ldag::Operator)::Vector{Float64}
    energy=[GaussianIntegral(Hermite(1, -dn[i+1], 0, i)*(Ldag*Lop*Hermite(1, -dn[i+1], 0, i))) for i=0:power]
    Float64.(nsolve.(SymPy.diff.(energy,dn),ones(power+1)))
end
calcu0(power::Integer, dn::Vector{SymPy.Sym}, Lop::Operator, Ldag::Nothing)::Vector{Float64}=ones(power+1);

function daef(out,du,u,Mclachlanfuncs::Vector{Function},t)
    for i=1:length(Mclachlanfuncs)
        out[i]::Float64=Mclachlanfuncs[i](du...,u...,t)::Float64
    end
end

function performMclachlan(power::Integer, Lop::Operator;tfinal::Float64,Ldag::Operator=nothing)::NamedTuple{(:Ks,:solution),Tuple{Vector{Sym},DAESolution{Float64}}}
    K::Vector{Sym}=[SymFunction("K$i")(t) for i in 0:power] #These are the coefficients
    #I need the widths as both functions of t and constants, so I can take their derivatives
    dt::Vector{Sym}=[SymFunction("dt$i")(t) for i in 0:power] #These are the widths/variational parameters as functions of t
    dn::Vector{Sym}=[symbols("dn$(i)") for i in 0:power]#These are the widths as constants
    ddt::Vector{Sym}=[symbols("ddt$i") for i in 0:power]

    println("Mathematica solving")
    psi=[Hermite(K[i+1], -dn[i+1], 0, i) for i=0:power]
    @time Ksol::Vector{Sym}=calcKsol(psi,power,K,Lop)

    #Now substitute in our solutions for the weights
    psisubs::Terms=Terms([subs(Ksol[i], dn, dt)*Hermite(1,-dt[i],0,power) for i=1:power+1])

    println("Substituted")
    #Setting up Mclachlan 2.5
    Theta::Terms=SymPy.diff(psisubs,t)
    LpsiMinusTheta::Terms=Lop*psisubs-Theta
    println("Integrating")
    @time Mclachlan::Vector{Sym}=Sym[SymPy.subs(SymPy.subs(GaussianIntegral(SymPy.diff(Theta,dt[i])*LpsiMinusTheta),SymPy.diff.(dt,Ref(t)),ddt),dt,dn) for i in 1:power+1]
    println("Lambdifying")
    Mclachlanfuncs::Vector{Function}=(eval ∘ lambdify).(Mclachlan,Ref(union(ddt,dn,[t])),invoke_latest=false,fns=Dict("Pow"=>:^))
    println("Solving for u0")
    @time u0=calcu0(power,dn,Lop,Ldag)

    #=possibledu0=solve(SymPy.subs.(SymPy.subs.(Mclachlansubs,Ref(Dict(dn[i]=>u0[i] for i in 1:length(dn)))),Ref(Dict(t=>0))),ddt)
    du0=[nothing for i=0:power]
    for i in possibledu0
        if !any(i.<=0)
            du0=[x for x in i]
            break
        end
    end=#
    println("Numerical solving")
    du0=zeros(power+1)
    tspan=(0.0,tfinal)
    f=DAEFunction(daef,syms=Symbol.(dn))
    prob=DAEProblem(f,du0,u0,tspan,Mclachlanfuncs,differential_vars=fill(true,power+1))
    return (Ks=Ksol,solution=@time solve(prob,daskr()))
end;
#For Fokker-Planck:

a, c, g = (5,1,2);

Kop=(a*X+c*X^3)+(g*Dx);
Lop=Dx*Kop;
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;

tfinal=15.0;

result=@time performMclachlan(0, Lop,tfinal=tfinal,Ldag=Ldag);
println(result)
(;Ks,solution)=result;
df=DataFrame(solution)

Kfunc=lambdify(Array(Ks),union([:t],Symbol.("dn$i" for i=0:length(Ks)-1)))
Knumerics=Kfunc.(solution.t,[solution[i,:] for i in 1:length(Ks)]...)

CSV.write(cwd*"/solution.csv",df)
CSV.write(cwd*"/Ks.csv",DataFrame(mapreduce(permutedims, vcat, Knumerics),:auto))

plot(solution)