using NLsolve
using SnoopCompile
using AbstractTrees
using ProfileView
using Profile

include("HermiteSymbolic.jl");
include("MathematicaSymbolics.jl");
include("JuliaMathematicaSymbolic.jl");

genPsi(power::Int, K::Vector{Num}, dn::Vector{Num},even::Bool=false)::Vector{Hermite}=begin
    if even
        return [Hermite(K[i+1], -dn[i+1], 0, 2i) for i=0:power]
    else
        return [Hermite(K[i+1], -dn[i+1], 0, i) for i=0:power]
    end
end

function calcKsol(Kpsi::Vector{Hermite},power::Int,K::Vector{Num},t::Num,Lop::Operator)::Vector{Num}
    #Put the L operator (it's a parameter of the function) into the basis of Hermites
    Lnm::Vector{Vector{Num}}=[[GaussianIntegral(Kpsi[n]*Lop(Kpsi[m]), Val(true)) for m in 1:power+1] for n in 1:power+1]
    #Mclachlan 2.4 left-hand side
    KDtK::Vector{Num}=[GaussianIntegral(derivative(Kpsi[i],t)*Kpsi[i],Val(true)) for i in 1:power+1]

    #Mclachlan 2.4
    Kequations::Vector{Equation}=Equation[KDtK[i]~sum(Lnm[i])[1] for i in 1:size(Lnm)[1]]

    #Solve them with initial condition of even weights
    #This is the step that freezes but that Mathematica can do
    DSolveMathematica(Kequations, K, t; ics=(substitute.(K,Ref(t=>0))::Vector{Num}).~fill(1/(power+1),length(K)))
end;

subsK(power::Int, Ksol::Vector{Num},dn::Vector{Num},even::Bool=false)::Terms=begin
    if even
        return Terms([Ksol[i]*Hermite(1,-dn[i],0,2(i-1)) for i=1:power+1])
    else
        return Terms([Ksol[i]*Hermite(1,-dn[i],0,(i-1)) for i=1:power+1])
    end
end

function genMclachlan(t::Num, dn::Vector{Num}, Lop::Operator, psisubs::Terms)
    Theta::Terms=substitute(derivative(psisubs,t),Dict(t=>0)) #Substitute t=0?
    LpsiMinusTheta::Terms=Lop*substitute(psisubs,Dict(t=>0))-Theta

    return @. GaussianIntegral(derivative($Ref(Theta),dn) * $Ref(LpsiMinusTheta),$Val{$length(Theta)*3+$length(LpsiMinusTheta)>500}())
end
function genMclachlan(t::Num, dn::Symbolics.Arr{Num, 1}, params::Vector{Num}, Lop::Operator, psisubs::Terms)
    Theta::Terms=@time "Theta" substitute(derivative(psisubs,t),Dict(t=>0)) #Substitute t=0?
    LpsiMinusTheta::Terms=@time "LpsiMinusTheta" Lop*substitute(psisubs,Dict(t=>100))-Theta

    retFuncs::Vector{Vector{Expr}} = fill(Vector{Expr}(undef, 3*length(Theta)*length(LpsiMinusTheta)), length(dn))
    for i=1:length(dn)
        retFuncs[i] .= GaussianIntegral(Function, derivative(Theta, dn[i]) * LpsiMinusTheta, dn, params...)
    end
    return retFuncs
end

#Change this to not have variables
function calcu0(power::Int, Lop::Operator, Ldag::Operator)::Vector{Float64}
    @variables localWidthsarr[1:power+1]
    localWidths::Vector{Num}=Symbolics.scalarize(localWidthsarr)
    energy::Vector{Num}=[GaussianIntegral(Hermite(1, -localWidths[i+1], 0, i)*(Ldag*Lop*Hermite(1, -localWidths[i+1], 0, i))) for i=0:power]
    mathematica_to_expr(weval(MathLink.WSymbol("Part")(MathLink.WSymbol("SolveValues")(expr_to_mathematica(@. (derivative(energy,localWidths)~0)),expr_to_mathematica(localWidths),MathLink.WSymbol("Reals")),1))::MathLink.WExpr)::Vector{Float64}
end
calcu0(power::Int, dn::Vector{Num}, Lop::Operator, Ldag::Nothing)::Vector{Float64}=ones(power+1);

function performMclachlan(power::Int, Lop::Operator, params::Vector{Num}, useEvens::Bool=false)#::NamedTuple{(:Ks,:expr),Tuple{Vector{Num},Vector{Num}}}
    #=@variables t::Float64
    K::Vector{Num}=Num[Symbolics.variable(:Karr,i,T=Symbolics.FnType{Tuple{Float64},Float64})(t) for i=1:power+1]
    dn::Vector{Num}=Num[Symbolics.variable(:dnarr,i,T=Float64) for i=1:power+1]
    dt::Vector{Num}=Num[Symbolics.variable(:dtarr,i,T=Symbolics.FnType{Tuple{Float64},Float64})(t) for i=1:power+1]=#
    
    @variables t Karr(t)[1:power+1] dtarr(t)[1:power+1] dnarr[1:power+1]
    K::Vector{Num}, dt::Vector{Num}, dn::Vector{Num} = Symbolics.scalarize.((Karr, dtarr, dnarr)::Tuple{Symbolics.Arr{Num, 1},Symbolics.Arr{Num, 1},Symbolics.Arr{Num, 1}})::NTuple{3, Vector{Num}}
    

    #println("Mathematica solving")
    psi::Vector{Hermite}=genPsi(power,K,dn,useEvens)

    #@time "Calculating K using Mathematica" Ksol::Vector{Num}=calcKsol(psi,power,K,t,Lop)
    Ksol::Vector{Num}=calcKsol(psi,power,K,t,Lop)
    #Now substitute in our solutions for the weights
    psisubs::Terms=subsK(power,Ksol,dn,useEvens)

    #println("Substituted")
    #Setting up Mclachlan 2.5
    #Mclachlan=@time "Mclachlan 2.5" genMclachlan(t, dn, Lop, psisubs)
    println("Generating Mclachlan functions")
    Mclachlan=genMclachlan(t, dn, Lop, psisubs)
    #subt0::Vector{Num} = @time substitute.(Mclachlan,Ref(t => 0))
    #withZeroDerivatives::Vector{Equation} = @time mathematica_to_expr(weval(MathLink.WSymbol("Replace")(expr_to_mathematica(Mclachlan .~ 0),MathLink.WSymbol("Rule")(MathLink.WSymbol("t"),0))))

    return (Ks=Ksol, expr=Mclachlan)
end;


#For Fokker-Planck:

params = @variables a c g;
#a, c, g = (5,1,3);

Kop=(a*X+c*X^3)+(g*Dx);
Lop=Dx*Kop;
Ldag=-(a*X+c*X^3)*Dx+g*Dx^2;

#=
(;Ks,expr) = performMclachlan(0, Lop, params, true);
tinf = expr
staleinstances(tinf)


tinf = @snoopi_deep performMclachlan(0, Lop, params, true);
ProfileView.view(flamegraph(tinf))
@profile performMclachlan(0, Lop, params, true)

import PyPlot
mref2, ax = pgdsgui(tinf)
=#

(;Ks,expr) = @time "Everything" performMclachlan(1, Lop, params, true);
result1 = performMclachlan(0, Lop, params, true);
mexpr1 = expr_to_mathematica(result1.expr);

@variables dnarr[1:1]
@time weval(W`FindRoot`(W`ReplaceAll`(mexpr1, expr_to_mathematica(Dict(a => 5, c => 1, g => 3))), W`List`(expr_to_mathematica(dnarr[1]), 1)))

mexpr = @time expr_to_mathematica(expr);
@variables dnarr[1:2]
dnsolve=@time mathematica_to_expr(weval(W`ReplaceAll`(expr_to_mathematica(Symbolics.scalarize(dnarr)),W`FindRoot`(W`ReplaceAll`(mexpr, expr_to_mathematica(Dict(a => 5, c => 1, g => 3))), W`List`(W`List`(expr_to_mathematica(dnarr[1]), 0.907),W`List`(expr_to_mathematica(dnarr[2]), 0.746))))))

Ksubs = substitute.(Ks,Ref(Dict(a => 5, c => 1, g => 3, Dict(Symbolics.scalarize(dnarr) .=> dnsolve)...)))
psisol = Hermite(Ksubs[1], -dnsolve[1], 0, 0)+Hermite(Ksubs[2], -dnsolve[2], 0, 2)
psisolnormed= psisol/GaussianIntegral(psisol)
substitute(psisolnormed, Dict(t=>100))
substitute(GaussianIntegral(psisolnormed*(substitute(Ldag*Lop*psisolnormed, Dict(a => 5, c => 1, g => 3)))),t=>145)


cfunc = build_function(expr,Symbolics.scalarize(dnarr),params,target=Symbolics.CTarget(), expression=Val(false));

open("C:/users/ethan/OneDrive/Desktop/julia_test/julia_func.cpp","w") do io
    write(io,cfunc)
end

funcs = Vector{Function}[]
function populate_funcs!(funcs::Vector{Vector{Function}}, expr::Vector{Vector{Expr}})
    for i=1:length(expr)
        i_funcs = Vector{Function}(undef, length(expr[i]))
        Threads.@threads for j=1:length(expr[i])
            i_funcs[j] = @time "Term $j" eval(expr[i][j])::Function
        end
        push!(funcs, i_funcs)
    end
    funcs
end
@time populate_funcs!(funcs, expr);

function call_vector_and_sum(funcs::Vector{Function}, args...)
    val1=funcs[1](args...)
    sums = zeros(typeof(val1), Threads.nthreads())
    sums[1]=val1
    Threads.@threads for i=2:length(funcs)
        sums[Threads.threadid()] += funcs[i](args...)
    end
    return sum(sums)
end
function call_vector_and_sum_Float64(funcs::Vector{Function}, args...)
    sums::Vector{Float64} = zeros(Float64, Threads.nthreads())
    Threads.@threads for i=1:length(funcs)
        sums[Threads.threadid()] += @time "Function $i" funcs[i](args...)
    end
    return sum(sums)
end
function fullFunc(theexpr::Vector{Vector{Function}}, args...)
    vec = zeros(length(theexpr))
    for i=1:length(theexpr)
        vec[i] = @time "dn $i" call_vector_and_sum_Float64(theexpr[i], args...)
    end
    vec
end

@time fullFunc(funcs, [1.0], 1.0, 1.0, 1.0)

#=

@variables dnarr[1:length(Ks)]

open("C:/Users/ethan/Downloads/Globus"*"/test_output.jl","w") do io
    write(io,"generated_expr_function(dnarr::Vector)=")
    write(io, string(expr)[4:end])
end;
include("C:/Users/ethan/Downloads/Globus"*"/test_output.jl");

param_vals= [5,1,3]

fn = @time eval(build_function(expr,Symbolics.scalarize(dnarr)...,params...,linenumbers=false,skipzeros=true)[1])
Cfn = @time build_function(expr,[dnarr;params],target=Symbolics.CTarget())
@report_opt eval(build_function(expr,Symbolics.scalarize(dnarr),linenumbers=false,skipzeros=true)[2])

@time nlsolve(fn, [1.3], autodiff=:forward)

@time weval(W`NSolveValues`(expr_to_mathematica(expr.~0),expr_to_mathematica(Symbolics.scalarize(dnarr)),W`Reals`))
@time (weval(W`NSolveValues`(expr_to_mathematica(substitute(expr.~0,Dict(zip(params,param_vals)))),expr_to_mathematica(Symbolics.scalarize(dnarr)),W`Reals`)))


calcu0(0, Lop, Ldag)
=#
#=
@named ns = NonlinearSystem(eqs, dnarr, params)
nlsys_func = eval(generate_function(ns)[2]) # second is the inplace version
j_func = eval(generate_jacobian(ns)[2]) # second is in-place


nlsolve((out,x)->nlsys_func(out,x,param_vals),(out,x)->j_func(out,x,param_vals),[1.3])

df=DataFrame(solution)

calcu0(length(Ks)-1, Symbolics.scalarize(dn), Lop, Ldag)

solveNonlinearSystem(ns, [5,1,2], ones(Float64, length))

Kfunc=eval(Base.remove_linenums!(build_function(Ks,[t;dn])[1]))
Knumerics=(Kfunc ∘ (x->[x...])).(collect(zip(solution.t,[solution[i,:] for i in 1:length(Ks)]...)))

CSV.write(cwd*"/solution.csv",df)
CSV.write(cwd*"/Ks.csv",DataFrame(mapreduce(permutedims, vcat, Knumerics),:auto))

plot(solution)=#