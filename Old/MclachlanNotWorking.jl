using DifferentialEquations
using DASKR
using Plots
using DataFrames
using CSV

const cwd=splitdir(Base.source_path())[1]
const daskrsolver=daskr();

include("HermiteNotWorking.jl");
include("JuliaMathematica.jl");

function MclachlanC(Kpsi::Vector{Hermite},power::Integer,K::Vector{SymPy.Sym},Lop::Operator)::Vector{SymPy.Sym}
    #Put the L operator (it's a parameter of the function) into the basis of Hermites
    Lnm::Vector{Vector{SymPy.Sym}}=[[GaussianIntegral(Kpsi[n]*Lop(Kpsi[m])) for m in 1:power+1] for n in 1:power+1]
    #Mclachlan 2.4 left-hand side
    KDtK::Vector{SymPy.Sym}=[GaussianIntegral(Dt(Kpsi[i])*Kpsi[i]) for i in 1:power+1]

    #Mclachlan 2.4
    Kequations::Vector{SymPy.Sym}=[simplify(SymPy.Eq(KDtK[i],sum(Lnm[i])[1])) for i in 1:size(Lnm)[1]]

    #Solve them with initial condition of even weights
    #This is the step that freezes but that Mathematica can do
    @time DSolveMathematica(Kequations, K; ics=Dict(K[i](0)=>1/sqrt(power+1) for i in 1:power+1))
end;

function calcu0(power::Integer, dn::Vector{SymPy.Sym}, Lop::Operator, Ldag::Operator)::Vector{Float64}
    energy::Vector{SymPy.Sym}=[GaussianIntegral(Hermite(1, -dn[i+1], 0, i)*(Ldag*Lop*Hermite(1, -dn[i+1], 0, i))) for i=0:power]
    println(energy)
    return energy
    dsol=solve(SymPy.diff.(energy,dn),Tuple(i for i in dn))
    u0=map(v->Float64(real(v[1])),filter.(Ref(x->abs(imag(x))<1e-3 && real(x)>0),dsol))
end;
calcu0(power::Integer, dn::Vector{SymPy.Sym}, Lop::Operator, Ldag::Nothing)::Vector{Float64}=ones(power+1);

function performMclachlan(power::Integer, Lop::Operator;tfinal=5,Ldag::Union{Nothing,Operator}=nothing)#::NamedTuple{(:Ks, :solution),Tuple{Vector{SymPy.Sym},DAESolution}}
    K::Vector{SymPy.Sym}=[SymFunction("K$i")(t) for i in 0:power] #These are the coefficients
    #I need the widths as both functions of t and constants, so I can take their derivatives
    dt::Vector{SymPy.Sym}=[SymFunction("dt$i")(t) for i in 0:power] #These are the widths/variational parameters as functions of t
    dn::Vector{SymPy.Sym}=[symbols("dn$(i)") for i in 0:power]#These are the widths as constants
    ddt::Vector{SymPy.Sym}=[symbols("ddt$i") for i in 0:power]
    
    #Switch all dn->dt or dt->dn
    dnTodt(herm) = subs(herm, dn, dt)
    dtTodn(herm) = subs(herm, dt, dn)

    #The Hermite basis state of a certain order
    Kpsi(power::Integer)=Hermite(K[power+1], -dn[power+1], 0, power)
    Kpsidt(power::Integer)=Hermite(K[power+1], -dt[power+1], 0, power)

    psi=Terms(sum(Kpsi,0:power))

    Ksol::Vector{SymPy.Sym}=MclachlanC(psi.terms,power,K,Lop)
    
    #Now substitute in our solutions for the weights
    psisubs::Terms=dnTodt(SymPy.subs(psi,Dict(K[i]=>Ksol[i] for i in 1:power+1)))

    println("Substituted")
    #Setting up Mclachlan 2.5
    Theta::Terms=Dt*psisubs
    LpsiMinusTheta::Terms=Lop*psisubs-Theta
    println("Integrating")
    @time Mclachlan::Vector{SymPy.Sym}=[GaussianIntegral(SymPy.diff(Theta,dt[i])*LpsiMinusTheta) for i in 1:power+1]
    @time Mclachlansubs::Vector{SymPy.Sym}=SymPy.subs.(SymPy.subs.(Mclachlan,Ref(SymPy.diff.(dt,Ref(t))),Ref(ddt)),Ref(dt),Ref(dn))
    Mclachlanfuncs=[lambdify(Mclachlansubs[i],union(ddt,dn,[t])) for i in 1:length(Mclachlansubs)]
    Mclachlanfunc(args...)::Vector{Float64}=begin
        [Float64(f(args...)) for f in Mclachlanfuncs] 
    end
    println("Setting up")
    function daef(out,du,u,p,t)
        result=Mclachlanfunc(du...,u...,t)
        for i in 1:length(result)
            out[i]=result[i]
        end
    end
    println("Function")
    
    u0=calcu0(power,dn,Lop,Ldag)
    return u0

    """possibledu0=solve(SymPy.subs.(SymPy.subs.(Mclachlansubs,Ref(Dict(dn[i]=>u0[i] for i in 1:length(dn)))),Ref(Dict(t=>0))),ddt)
    du0=[nothing for i=0:power]
    for i in possibledu0
        if !any(i.<=0)
            du0=[x for x in i]
            break
        end
    end"""
    du0=zeros(power+1)
    tspan=(0.0,tfinal)
    f=DAEFunction(daef,syms=Symbol.(dn))
    prob=DAEProblem(f,du0,u0,tspan,differential_vars=fill(true,power+1))
    @time solution=DifferentialEquations.solve(prob,daskrsolver)::DAESolution
    return (Ks=Ksol::Vector{SymPy.Sym},solution=solution::DAESolution)
end;
#For Fokker-Planck:

a, c, g = (5,1,2);

Kop=(a*X+c*X^3)+(g*Dx);
Lop=Dx*Kop;
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;

tfinal=15;

result=@time performMclachlan(0, Lop,tfinal=tfinal,Ldag=Ldag);
@code_warntype performMclachlan(0, Lop,tfinal=tfinal,Ldag=Ldag)
(;Ks,solution)=result;
df=DataFrame(solution)

Kfunc=lambdify(Array(Ks),union([:t],Symbol.("dn$i" for i=0:length(Ks)-1)))
Knumerics=Kfunc.(solution.t,[solution[i,:] for i in 1:length(Ks)]...)

CSV.write(cwd*"/solution.csv",df)
CSV.write(cwd*"/Ks.csv",DataFrame(mapreduce(permutedims, vcat, Knumerics),:auto))

tsol=solution(tfinal)
tK=Kfunc(t,tsol...).subs(Dict(t=>tfinal))
tPsi=sum(Hermite(tK[i],-tsol[i],0,i-1) for i=1:length(tsol))
tPsinorm=tPsi/GaussianIntegral(tPsi)
plot(Sym(tPsinorm))