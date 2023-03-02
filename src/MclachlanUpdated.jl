include("HermiteSymbolic.jl");
include("../JuliaMathematica/SymbolicsMathLink.jl");

genPsi(power::Int, dn::Vector{Num},even::Bool=false)= begin
    if even
        return Hermite[Hermite(1, -dn[i+1], 0, 2i) for i=0:power]
    else
        return Hermite[Hermite(1, -dn[i+1], 0, i) for i=0:power]
    end
end

function calcKsol(MclachlanX::Terms, psi::Vector{Hermite},power::Integer,K::Vector{Num},t::Num,Lop::Operator)

    #Mclachlan 2.4
    Kequations::Vector{Equation}=[GaussianIntegral(MclachlanX*psi[i])~0 for i=1:power+1]

    #Solve them with initial condition of even weights
    #This is the step that freezes but that Mathematica can do
    initconds = [1;zeros(Int,length(K)-1)]#fill((1/(power+1))::Float64,length(K))#
    DSolveMathematica(Kequations, K, t; ics=(substitute.(K,Ref(t=>0))::Vector{Num}).~initconds)
end;

subsK(power::Int, Ksol::Vector{Num},dn::Vector{Num},even::Bool=false)::Terms=begin
    if even
        return Terms([Ksol[i]*Hermite(1,-dn[i],0,2(i-1)) for i=1:power+1])
    else
        return Terms([Ksol[i]*Hermite(1,-dn[i],0,(i-1)) for i=1:power+1])
    end
end;

function genMclachlan(dn::Vector{Num}, Theta::Terms, MclachlanX::Terms)
    return @. GaussianIntegral(derivative($Ref(Theta),dn) * $Ref(MclachlanX),$Val(true))
end;

function performMclachlan(power::Int, Lop::Operator, params::Vector{Num}, useEvens::Bool=false)
    @variables t Karr(t)[1:power+1] dnarr[1:power+1]
    K::Vector{Num}, dn::Vector{Num}, = scalarize.((Karr, dnarr)::NTuple{2,Symbolics.Arr{Num, 1}})::NTuple{2, Vector{Num}}

    psi::Vector{Hermite}=@time "genPsi" genPsi(power,dn,useEvens)

    Theta::Terms=@time "Theta" Terms(@. derivative(K, t)*psi)
    MclachlanX::Terms = @time "MclachlanX" Theta-sum(@. Lop(K*psi))
    
    Ksol::Vector{Num}=@time "calcKsol" calcKsol(MclachlanX, psi,power,K,t,Lop)
    #Now substitute in our solutions for the weights
    #psisubs::Terms=subsK(power,Ksol,dn,useEvens)

    #println("Substituted")
    #Setting up Mclachlan 2.5
    substitution = Dict(K.=>Ksol)
    Mclachlan=@time "genMclachlan" genMclachlan(dn,(@time "substituteTheta" substitute(Theta, substitution)),(@time "substituteMclachlanX" substitute(MclachlanX,substitution)))

    return (Ks=Ksol, expr=Mclachlan, dn=dn)
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

result=@time performMclachlan(1, Lop, params, true);

mdn = expr_to_mathematica(result.dn)

mexpr = @time W"Thread"(W"Equal"(expr_to_mathematica(result.expr),0));
mexprWithParams=@time weval(W"ReplaceAll"(mexpr,W`t->10`); Pi=pi, Dict(Symbol.(params).=>param_vals)...);
mexprWithParams
solution = weval(W"NSolve"(mexprWithParams,mdn,W"Reals"))