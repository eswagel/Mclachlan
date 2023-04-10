using Symbolics
using MathLink

const Mtypes = Union{MathLink.WTypes,Int8,Int16,Int32,Int64,Int128,UInt8,UInt16,UInt32,UInt64,UInt128,Float16,Float32,Float64,ComplexF16,ComplexF32,ComplexF64,Rational}

function decode_piecewise(lists::Vector{Vector{Num}}, lastval)
    @nospecialize
    second_to_lastval::Vector = lists[end]
    ret = Symbolics.IfElse.ifelse(second_to_lastval[2], second_to_lastval[1], lastval)
    for i=length(lists)-2:-1:1
        cur_val::Vector = lists[i]
        ret = Symbolics.IfElse.ifelse(curval[2], curval[1], ret)
    end
    ret
end
function decode_piecewise(lists::Vector{Num}, lastval)
    @nospecialize
    Symbolics.IfElse.ifelse(lists[2], lists[1], lastval)
end
function decode_piecewise(lists::Vector{Vector})
    @nospecialize
    ret = Symbolics.IfElse.ifelse(lists[1][2], lists[1][1], nothing)
    for i=2:length(lists)
        ret = Symbolics.IfElse.ifelse(lists[i][2], lists[i][1], ret)
    end
    ret
end

numize_if_not_vector(x::Vector{Num})=x
numize_if_not_vector(x::Vector)=numize_if_not_vector.(x)
numize_if_not_vector(x::Tuple)=[x...]
numize_if_not_vector(x::Num)=x
numize_if_not_vector(x)=Num(x)
numize_if_not_vector(x::Pair)=numize_if_not_vector(x[1]) => numize_if_not_vector(x[2])

Num(x::Complex)=Num(real(x))+Num(imag(x))*im


# Define a dictionary that maps Mathematica function names to their Julia equivalents
const MATHEMATICA_TO_JULIA_FUNCTIONS::Dict{String,Function} = Dict(
    "Plus" => +,
    "Times" => *,
    "Power" => ^,
    "Subtract" => -,
    "Divide" => /,
    "Sin" => sin,
    "Cos" => cos,
    "Tan" => tan,
    "ArcSin" => asin,
    "ArcCos" => acos,
    "ArcTan" => atan,
    "Sinh" => sinh,
    "Cosh" => cosh,
    "Tanh" => tanh,
    "ArcSinh" => asinh,
    "ArcCosh" => acosh,
    "ArcTanh" => atanh,
    "Exp" => exp,
    "Log" => log,
    "Sqrt" => sqrt,
    "Abs" => abs,
    "Rational" => //,
    "List" => function (x...)
        [numize_if_not_vector.(x)...]
    end,
    "Rule" =>(a,b)->a=>b,
    "C" => i->Symbolics.variable("C_$i"),
    "Part" => getindex,
    "Equal" => Base.:~,
    "Piecewise" => decode_piecewise,
    "Greater" => >,
    "Less" => <,
    "GreaterEqual" => >=,
    "LessEqual" => <=,
    "Complex" => (a,b)->a+b*im,
)

const JULIA_FUNCTIONS_TO_MATHEMATICA = Dict(
    :+ => "Plus",
    :* => "Times",
    :^ => "Power",
    :- => "Subtract",
    :/ => "Divide",
    :sin => "Sin",
    :cos => "Cos",
    :tan => "Tan",
    :asin => "ArcSin",
    :acos => "ArcCos",
    :atan => "ArcTan",
    :sinh => "Sinh",
    :cosh => "Cosh",
    :tanh => "Tanh",
    :asinh => "ArcSinh",
    :acosh => "ArcCosh",
    :atanh => "ArcTanh",
    :exp => "Exp",
    :log => "Log",
    :sqrt => "Sqrt",
    :abs => "Abs",
    :// => "Rational",
    :~ => "Equal",
    :> => "Greater",
    :< => "Less",
    :>= => "GreaterEqual",
    :<= => "LessEqual",
    :!= => "Unequal",
)

expr_to_mathematica(expr::Expr)::Mtypes=expr_to_mathematica(expr.head,expr.args)
#=expr_to_mathematica(function_name::String,args::Vector)=begin
    println(function_name,": ",args)
    if haskey(JULIA_FUNCTIONS_TO_MATHEMATICA, Symbol(function_name))
        return W"$(JULIA_FUNCTIONS_TO_MATHEMATICA[Symbol(function_name)])"(expr_to_mathematica.(args)...)
    else
        return W"$(function_name)"(expr_to_mathematica.(args)...)
    end
end=#
expr_to_mathematica_differential_checker(function_head::Symbol,args::Vector,::Nothing)::MathLink.WExpr=W"$(string(function_head))"(expr_to_mathematica.(args)...)
expr_to_mathematica_differential_checker(function_head::Symbol,args::Vector,m::RegexMatch)::MathLink.WExpr=begin
    return MathLink.WSymbol("D")(expr_to_mathematica.(args)...,MathLink.WSymbol("$(m[1])"))
end
expr_to_mathematica_vector_handler(expr::MathLink.WSymbol,index::Integer)=MathLink.WSymbol("$(expr.name)■$(index)")
expr_to_mathematica_vector_handler(expr::MathLink.WExpr,index::Integer)=begin
    MathLink.WSymbol("$(expr.head.name)■$(index)")(expr.args...)
end
expr_to_mathematica(function_head::Symbol,args::Vector)::Mtypes=begin
    if function_head==:call
        return expr_to_mathematica(Symbol(args[1]),args[2:end])
    elseif function_head==:inv
        return MathLink.WSymbol("Power")(expr_to_mathematica(args[1]),-1)
    elseif function_head==:getindex
        would_be = expr_to_mathematica(args[1])
        return expr_to_mathematica_vector_handler(would_be, args[2])
    elseif function_head==:if
        cond = expr_to_mathematica(args[1])
        ifval = expr_to_mathematica(args[2])
        elseval = expr_to_mathematica(args[3])
        return MathLink.WSymbol("Piecewise")(MathLink.WSymbol("List")(MathLink.WSymbol("List")(ifval, cond)),elseval)
    elseif haskey(JULIA_FUNCTIONS_TO_MATHEMATICA, function_head)
        return MathLink.WSymbol("$(JULIA_FUNCTIONS_TO_MATHEMATICA[function_head])")(expr_to_mathematica.(args)...)
    else
        fstring=string(function_head)
        m = match(r"Differential\(([^)]*)\)",fstring)
        return expr_to_mathematica_differential_checker(function_head,args,m)
    end
end
expr_to_mathematica_symbol_vector_checker(str::String,m::Nothing)=MathLink.WSymbol(str)
expr_to_mathematica_symbol_vector_checker(str::String,m::RegexMatch)=begin
    replacements=("₁"=>"1","₂"=>"2","₃"=>"3","₄"=>"4","₅"=>"5","₆"=>"6","₇"=>"7","₈"=>"8","₉"=>"9","₀"=>"0")
    MathLink.WSymbol("$(m[1])■$(replace(m[2],replacements...))")
end
expr_to_mathematica(symbol::Symbol)::MathLink.WSymbol=begin
    if haskey(JULIA_FUNCTIONS_TO_MATHEMATICA, symbol)
        return MathLink.WSymbol(JULIA_FUNCTIONS_TO_MATHEMATICA[symbol])
    else
        symString=string(symbol)
        m=match(r"([^₁₂₃₄₅₆₇₈₉₀]*)([₁|₂|₃|₄|₅|₆|₇|₈|₉|₀]+)$",symString)
        return expr_to_mathematica_symbol_vector_checker(symString,m)
    end
end
(expr_to_mathematica(num::T)::T) where {T<:Mtypes}=num
expr_to_mathematica(eq::Equation)::MathLink.WExpr=MathLink.WSymbol("Equal")(expr_to_mathematica(Symbolics.toexpr(eq.lhs)::Union{Expr, Symbol, Int, Float64, Rational}),expr_to_mathematica(Symbolics.toexpr(eq.rhs)::Union{Expr, Symbol, Int, Float64, Rational}))
(expr_to_mathematica(vect::Vector{T})::MathLink.WExpr) where T=MathLink.WSymbol("List")(expr_to_mathematica.(vect)...)
(expr_to_mathematica(mat::Matrix{T})::MathLink.WExpr) where T = expr_to_mathematica([mat[:,i] for i in 1:size(mat,2)])
expr_to_mathematica(num::Num)::Mtypes=expr_to_mathematica(Symbolics.toexpr(num)::Union{Expr, Symbol, Int, Float64, Rational})
expr_to_mathematica(dict::Dict)::MathLink.WExpr=begin
    rules = MathLink.WExpr[]
    for (key, val) in dict
        push!(rules, MathLink.WSymbol("Rule")(expr_to_mathematica(key), expr_to_mathematica(val)))
    end
    MathLink.WSymbol("List")(rules...)
end
expr_to_mathematica(sym::Symbolics.Symbolic)::MathLink.WExpr=begin
    expr::Expr = Symbolics.toexpr(sym)::Expr
    expr_to_mathematica(expr)
end
expr_to_mathematica(st::AbstractString)::MathLink.WSymbol=MathLink.WSymbol(st)
expr_to_mathematica(x::BigFloat)=Float64(x)
expr_to_mathematica(x::Irrational)=Float64(x)

expr_to_mathematica(num::Num, symbolic::Val{true})::MathLink.WExpr=begin
    expr::Expr = Symbolics.toexpr(num)::Expr
    expr_to_mathematica(expr,Val(true))
end
expr_to_mathematica(vect::Vector{Num}, symbolic::Val{true})::MathLink.WExpr=MathLink.WSymbol("List")(expr_to_mathematica.(vect, Val(true))...)
expr_to_mathematica(expr::Expr, symbolic::Val{true})::MathLink.WExpr=expr_to_mathematica(expr.head,expr.args)

mathematica_to_expr_vector_checker(head::AbstractString,args::Vector,::Nothing)=begin
    varname=Symbol(head)
    (vars,)=Symbolics.@variables varname
    vars
end
mathematica_to_expr_vector_checker(head::AbstractString,args::Vector,m::RegexMatch)=begin
    varname=Symbol(m[1])
    (vars,)=scalarize.(Symbolics.@variables $varname(mathematica_to_expr.(args)...)[1:parse(Int8,m[2])])
    vars[parse(Int8,m[2])]
end
mathematica_to_expr_differential_checker(head::MathLink.WExpr,args::Vector)::Num=begin
    if head.head.head.name=="Derivative"
        return (Differential(Symbolics.variable(args[1]))^head.head.args[1])(eval(Symbol(head.args[1])))
    else
        throw(ArgumentError("Not a derivative"))
    end
end
mathematica_to_expr_differential_checker(head::MathLink.WSymbol,args::Vector)=begin
    if head.name=="Power" && isa(args[1],MathLink.WSymbol) && args[1].name=="E"
        return exp(mathematica_to_expr.(args[2:end])...)
    elseif haskey(MATHEMATICA_TO_JULIA_FUNCTIONS, head.name)
        return MATHEMATICA_TO_JULIA_FUNCTIONS[head.name](mathematica_to_expr.(args)...)
    else
        m=match(r"(.+)■([0-9]+)",head.name)
        return mathematica_to_expr_vector_checker(head.name,args,m)
        #return mathematica.head.name,getproperty.(mathematica.args,Ref(:name))...)
    end
end
mathematica_to_expr_differential_checker(head,args::Tuple)=mathematica_to_expr_differential_checker(head,[args...])
mathematica_to_expr(mathematica::MathLink.WExpr)=mathematica_to_expr_differential_checker(mathematica.head,mathematica.args)
mathematica_to_expr(symbol::MathLink.WSymbol)=mathematica_to_expr(symbol,match(r"(.+)■([0-9]+)",symbol.name))
mathematica_to_expr(symbol::MathLink.WSymbol,::Nothing)=Symbolics.value(Symbolics.variable(symbol.name))
mathematica_to_expr(symbol::MathLink.WSymbol,m::RegexMatch)=begin
    varname=Symbol(m[1])
    (vars,)=scalarize.(Symbolics.@variables $varname[1:parse(Int8,m[2])])
    vars[parse(Int8,m[2])]
end
mathematica_to_expr(mathematica::MathLink.WReal)=mathematica_to_expr(weval(W"N"(mathematica)))
mathematica_to_expr(num::T) where T<:Number=begin
    rounded = round(num)
    if abs(rounded-num)< 1e-10
        return rounded
    else
        return num
    end
end

include("DSolveMathematica.jl")

function wcall(head::AbstractString, args...; kwargs...)
    return wcall(head, expr_to_mathematica.(args)...;kwargs...)
end
wcall(head::AbstractString, args::Vararg{Mtypes}; kwargs...) = mathematica_to_expr(weval(MathLink.WSymbol(head)(args...; kwargs...)))