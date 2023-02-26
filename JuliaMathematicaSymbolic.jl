using MathLink

DSolveMathematica(expr::AbstractVector{Num}, vars, args...; kwargs...)::Vector{Num}=DSolveMathematica(expr.~0, vars, args...; kwargs...)
DSolveMathematica(expr, vars::Num, args...; kwargs...)::Vector{Num}=DSolveMathematica(expr, [vars], args...; kwargs...)
DSolveMathematica(expr,vars::Symbolics.Arr{Num, 1}, args...;kwargs...)::Vector{Num}=DSolveMathematica(expr,[vars[i] for i=1:length(vars)],args...;kwargs...)
DSolveMathematica(expr::Equation,vars::AbstractVector{Num},args...; kwargs...)::Vector{Num}=DSolveMathematica([expr],vars,args...;kwargs...)
DSolveMathematica(expr::AbstractVector{Equation},vars::AbstractVector{Num}, t::Num; kwargs...)::Vector{Num}=DSolveMathematica(expr,vars,string(t);kwargs...)
DSolveMathematica(expr::AbstractVector{Equation}, vars::AbstractVector{Num}, t::String; ics::AbstractVector{Equation}=Equation[])::Vector{Num}=begin
    mexpr=expr_to_mathematica(expr)
    exprics=expr_to_mathematica(ics)
    mathExpr::MathLink.WExpr=(MathLink.WSymbol("DSolveValue")(weval(MathLink.WSymbol("Join")(mexpr,exprics)),expr_to_mathematica(vars),MathLink.WSymbol("$t")))
    dsolved::MathLink.WExpr=weval(mathExpr)
    mathematica_to_expr(dsolved)
end


