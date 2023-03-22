include("HermiteSymbolic.jl");
include("../JuliaMathematica/SymbolicsMathLink.jl");
using DifferentialEquations
import DASKR: daskr
#import Sundials: IDA
using Plots

genPsi(power::Int, dn::Vector{Num},skip::Int=1)=begin
    ret = Vector{Terms}(undef,power+1)
    for i=0:power
        ret[i+1] = Terms(Hermite(1,-dn[i+1],0,i*skip))#hermite_poly(skip*i,X)*Hermite(1, -dn[i+1], 0, 0)
    end
    ret
end

psi=genPsi(1,[dn[i] for i=1:length(dn)],2)
psinorm=[psi[i]/sqrt(GaussianIntegral(psi[i]^2)) for i=1:length(psi)]

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

function calcu0(H::Operator)
    @variables dn
    mdn = W"dn"
    psi = Hermite(1, -dn, 0, 0)

    energy = GaussianIntegral(psi*(H*psi))

    wcall("SolveValues",Differential(dn)(energy)~0,dn,"Reals")
end;

params = @variables a c g;

param_vals = [5.0, 1.0, 1.0]
#a, c, g = (5,1,2);

Lop=Dx*(a*X+c*X^3)+(g*Dx^2);
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;

tfinal=10.0;

power,skip = 1,2
result=@time performMclachlan(power, Lop, params, skip);

resultfunc=eval.(build_function(substitute(result.expr,Dict(Differential(t).([result.K;result.dt]).=>[result.dK;result.ddt])),[result.dK;result.ddt],[result.K;result.dt],params,t))

initcondsguess=[[[result.K[1],1.0]];[[result.K[i],1.0] for i=2:power+1];[[result.dt[i],i>1 ? 3.0 : 3.0] for i=1:power+1]]
initconds=substitute([result.K;result.dt],Dict(wcall("FindRoot",resultfunc[1](fill(-0.01,2(power+1)),[result.K;result.dt],param_vals,0.0),initcondsguess)))

daef=DAEFunction(resultfunc[2],syms=[:K,:d],paramsyms=[:a,:c,:g])
daep=DAEProblem(daef,fill(-0.01,2(power+1)),Symbolics.value.(initconds),[0.0,tfinal],param_vals)
daesol=solve(daep,daskr());

ksol=daesol[end][1:power+1]
dsol=daesol[end][power+2:end]

#hermitesol=sum([hermite_poly(2i,X)*Hermite(ksol[i+1],-dsol[i+1],0,0) for i=0:power]);
hermitesol=sum([Hermite(ksol[i+1],-dsol[i+1],0,skip*i) for i=0:power]);
hermitesolnorm=hermitesol/sqrt(GaussianIntegral(hermitesol*hermitesol));

substitute(GaussianIntegral(hermitesolnorm*(Ldag*Lop*hermitesolnorm)),Dict(params.=>param_vals))

plot(eval(build_function(convert(Num,hermitesolnorm),x)),-5,5)

#=
mdt = expr_to_mathematica(result.dt)
mK = expr_to_mathematica(result.K)

mexpr = @time W"Thread"(W"Equal"(expr_to_mathematica(result.expr),0));
mexprWithParams=@time weval(mexpr; Pi=pi, Dict(Symbol.(params).=>param_vals)...);

initconds = weval(expr_to_mathematica([result.K.~[1;zeros(Int,length(result.K)-1)];result.dt.~[3.0,3.0]]),t=0)

#solution = weval(W"NSolve"(mexprWithParams,mdn,W"Reals"))
msol=weval(W"NDSolveValue"(W"Join"(mexprWithParams,initconds),W"Join"(mdt,mK),W"List"(W"t",0,100),W`Method->{"EquationSimplification"->"Residual"}`));

msol1=weval(W"Part"(msol,1));
weval(`$msol1[2]`)
=#