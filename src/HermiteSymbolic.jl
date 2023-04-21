# Help me comment the code below

using SpecialFunctions
using Latexify
import Latexify: latexify
using Symbolics
import Symbolics: substitute, derivative#So I can overload them
using MathLink

#Built-in iszero is bad for Num
iszero(x::Num)=isprimitivetype(typeof(x.val)) && x.val==0
scalarize(arr::Symbolics.Arr{Num,1})::Vector{Num}=Num[arr[i] for i=1:length(arr)]
Symbolics.derivative(::Union{Int8, Int16, Int32, Int64, Int128, UInt8, UInt16, UInt32, UInt64, UInt128, Float16, Float32, Float64, ComplexF16, ComplexF32, ComplexF64},args...)=0

#Abstract type that will include Hermite and Terms
abstract type HermiteType end

#This struct models a Hermite basis state
struct Hermite <: HermiteType
    factor::Num #Coefficient in front
    a::Num #Coefficient of x^2 in the exponent (should be negative)
    b::Num #Coefficient of x in the exponent
    order::UInt8 #Power of x in the coefficient
end

#Single Hermite basis state
Hermite(a::Num,b::Num,order::Int)=Hermite(1,a,b,order)

#When substituting into a Hermite, substitute into each part 
Symbolics.substitute(herm::Hermite, args...)::Hermite=Hermite(Symbolics.substitute(herm.factor, args...), Symbolics.substitute(herm.a, args...), Symbolics.substitute(herm.b, args...), herm.order)

#From two arrays, create a Dict for substitution and substitute
Symbolics.substitute(herm::Hermite, syms1::Union{Vector{Num},Tuple{Num}}, syms2::Union{Vector{Num},Tuple{Num}})::Hermite = substitute(herm, Dict(syms1.=>syms2))
Symbolics.substitute(herm::Num, syms1::Union{Vector{Num},Tuple{Num}}, syms2::Union{Vector{Num},Tuple{Num}})::Num = substitute(herm, Dict(syms1.=>syms2))

Base.convert(::Type{Num},herm::Hermite)=herm.factor*x^herm.order*exp(herm.a*x^2+herm.b*x)

#Take the derivative of a single Hermite with respect to some symbol
#Always returns a Terms with three Hermites
Symbolics.derivative(herm::Hermite, var::Num)::Terms=begin
    result=Vector{Hermite}(undef,3) #Eventual output

    facdiff=derivative(herm.factor,var)
    
    result[1]=Hermite(facdiff,herm.a,herm.b,herm.order)

    adiff=derivative(herm.a,var)
    result[2]=Hermite(herm.factor*adiff,herm.a,herm.b,herm.order+2)

    bdiff=derivative(herm.b,var)
    result[3]=Hermite(herm.factor*bdiff,herm.a,herm.b,herm.order+1)

    return Terms(result)
end

#Struct that holds an array of Hermites
struct Terms <: HermiteType
    terms::Array{Hermite,1}
end

Terms(herm::Hermite)=Terms([herm])
Terms(terms::Terms)=terms
Base.zero(::Type{Terms})=Terms(Vector{Hermite}(undef,0))
Base.zero(::Terms)=Terms(Vector{Hermite}(undef,0)) #Zero is an empty array
Base.zero(::Hermite)=Hermite(0,0,0,0)
Base.zero(::Type{Hermite})=Hermite(0,0,0,0)

#Not working but would be convenient to be able to broadcast over terms
Base.length(terms::Terms)=length(terms.terms)
Broadcast.broadcast(func, terms::Terms) = Terms(Broadcast.broadcast(func, terms.terms))

#Derivative of terms = sum of derivatives
Symbolics.derivative(terms::Terms, args...)::Terms=begin
    sum=Vector{Hermite}(undef,length(terms)*3)
    for i=1:length(terms)
        @inbounds sum[(i-1)*3+1:i*3].=derivative(terms.terms[i], args...).terms
    end
    Terms(sum)
end

Symbolics.substitute(terms::Terms,args...)::Terms=begin
    result=Vector{Hermite}(undef,length(terms))
    for i=1:length(terms)
        @inbounds result[i]=substitute(terms.terms[i],args...)
    end
    return Terms(result)
end

Base.convert(::Type{Num},terms::Terms)=sum(convert.(Ref(Num),terms.terms))
Latexify.latexify(herm::HermiteType,args...;kwargs...)=Latexify.latexify(convert(Num,herm),args...;kwargs...)

#######Group algebra
#Checks if two Hermite basis states can be added together to be one, returns a Bool
addable(herm1::Hermite,herm2::Hermite)::Bool=isequal(herm1.a,herm2.a) && isequal(herm1.b,herm2.b) && isequal(herm1.factor,herm2.factor) && isequal(herm1.order,herm2.order)

#Multiplication
Base.:*(herm1::Hermite,herm2::Hermite)::Hermite=Hermite(herm1.factor*herm2.factor,herm1.a+herm2.a,herm1.b+herm2.b,herm1.order+herm2.order)
Base.:*(factor, herm2::Hermite)::Hermite=Hermite(factor*herm2.factor,herm2.a,herm2.b,herm2.order)
Base.:*(herm1::Hermite,factor)::Hermite=factor*herm1

#Division
Base.:/(herm::Hermite, fact)=herm*(1/fact)

#Exponentiation
Base.:^(herm::Hermite,p)::Hermite=Hermite(herm.factor^p,herm.a*p,herm.b*p,herm.order^p)

#Addition
Base.:+(herm1::Hermite,herm2::Hermite)::Terms=begin
    if iszero(herm1.factor) && iszero(herm2.factor)
        return zero(Terms)
    elseif addable(herm1,herm2)
        return Terms(Hermite(herm1.factor+herm2.factor, herm1.a, herm1.b, herm1.order))
    elseif iszero(herm1.factor)
        return Terms(herm2)
    elseif iszero(herm2.factor)
        return Terms(herm1)
    else
        return Terms([herm1,herm2])
    end
end
Base.:+(herm::Hermite,terms::Terms)::Terms=begin
    if iszero(herm.factor)
        return Terms(terms.terms)
    end
    return Terms([terms.terms;herm])
end
Base.:+(terms::Terms,herm::Hermite)=herm+terms

Base.:+(a::Terms,b::Terms)=begin
    return Terms([a.terms;b.terms])
end

#Subtraction

Base.:-(a::Hermite,b::Terms)::Terms=begin
    sum=Vector{Hermite}(undef,length(b)+1)
    sum[1]=a
    @inbounds @. sum[2:end]=-1*b.terms
    return sum
end
Base.:-(a::Terms,b::Hermite)::Terms=begin
    sum=Vector{Hermite}(undef,length(a)+1)
    @inbounds sum[1:end-1].=a.terms
    @inbounds sum[end]=-b
    return Terms(sum)
end
Base.:-(a::Terms,b::Terms)::Terms=begin
    sum=Vector{Hermite}(undef,length(a)+length(b))
    @inbounds sum[1:length(a)].=a.terms
    @inbounds sum[length(a)+1:end].=-1 .*b.terms
    return Terms(sum)
end
Base.:-(a::HermiteType)=(-1)*a


#Multiplying with Terms is just distribution
#Multiplication is multithreaded
Base.:*(a::Terms,b::Terms)::Terms=begin
    new_terms=Vector{Hermite}(undef,length(a)*length(b))
    L=length(b)
    Threads.@threads for i=1:length(a)
        for j=1:L
            @inbounds new_terms[(i-1)*L+j]=a.terms[i]*b.terms[j]
        end
    end
    return Terms(new_terms)
end
Base.:*(terms::Terms,herm::Hermite)::Terms=begin
    new_terms=Vector{Hermite}(undef,length(terms))
    for i=1:length(terms)
        @inbounds new_terms[i]=terms.terms[i]*herm
    end
    Terms(new_terms)
end
Base.:*(herm::Hermite,terms::Terms)=terms*herm
Base.:*(fact,terms::Terms)=Terms(fact .* terms.terms)
Base.:*(terms::Terms,fact)=fact*terms

#Division
Base.:/(terms::Terms, fact)=terms*(1/fact)

#Exponentiation
Base.:^(terms::Terms, p::Integer)=begin
    new_t = terms
    for i=1:p-1
        new_t = new_t*terms
    end
    new_t
end

function HalfIntegerGamma(n::Int)
    """Gamma function of n+1/2"""
    Float64(factorial(big(2n))/(4^n * factorial(big(n)))) * sqrt(pi)
end

#Integration of a Hermite basis state with a, b, and order.
#Can only have order up to 6 if b is not zero
GaussianIntegral(herm::Hermite)::Num=begin
    if iszero(herm.factor)
        return zero(Num)
    else
        a = herm.a
        b = herm.b
        order = herm.order
        if iszero(b)
            if isodd(order)
                return zero(Num)
            else
                return herm.factor*(-a)^(-(1//2) - order//2)*HalfIntegerGamma(Int(order/2))
            end
        elseif order==0
            result=1/sqrt(-a)
        elseif order==1
            result=-(b/(2*sqrt(-a)*a))
        elseif order==2
            result=-((1 - b^2/(2a))/(2*sqrt(-a)*a))
        elseif order==3
            result=-((3b*(1 - b^2/(6a)))/(4*(-a)^(3/2)*a))
        elseif order==4
            result=(3*(1 - b^2/a + b^4/(12a^2)))/(4*sqrt(-a)*a^2)
        elseif order==5
            result=-((15b*(1 - b^2/(3a) + b^4/(60a^2)))/(8*(-a)^(5/2)*a))
        elseif order==6
            result=-((15*(1 - (3b^2)/(2a) + b^4/(4a^2) - b^6/(120a^3)))/(8*sqrt(-a)*a^3))
        else
            throw(DomainError(order,"Integral for this order Hermite polynomial has not been calculated."))
        end
        return herm.factor*result*sqrt(pi)*exp(-b^2/(4a))
    end
end

#Integration of a Terms object (parallelized)
GaussianIntegral(terms::Terms)::Num=begin
    ints=zeros(Num,Threads.nthreads())
    Threads.@threads for i=1:length(terms.terms)
        @inbounds ints[Threads.threadid()]+=GaussianIntegral(terms.terms[i])
    end
    return sum(ints)
end
GaussianIntegral(terms::Terms,::Val{false})::Num=GaussianIntegral(terms)
#To not parallelize, use Val(true)
GaussianIntegral(terms::Terms,::Val{true})::Num=begin
    ints=Num(0)
    for i=1:length(terms.terms)
        @inbounds ints+=GaussianIntegral(terms.terms[i])
    end
    return ints
end

#An Operator manipulates the Hermites
struct Operator
    func::Function #Must be a function that takes a Hermite and returns a Terms
end

#Multiplying an operator by a Hermite applies it
Base.:*(op::Operator, herm::Hermite)::Terms=op.func(herm)::Terms

#Choose whether to parallelize based on the number of terms in the Terms object
operator_multiplication_parallel_handler(op::Operator, terms::Terms, ::Val{true})::Terms=begin
    sum=fill(Hermite[],Threads.nthreads())
    for i=1:length(terms.terms)
        sum[Threads.threadid()]=vcat(sum[Threads.threadid()],op.func(terms.terms[i]).terms)
    end
    return Terms(collect(Iterators.flatten(sum)))
end
operator_multiplication_parallel_handler(op::Operator, terms::Terms,::Val{false})::Terms=begin
    sum=zero(Terms)
    for i=1:length(terms.terms)
        sum+=op.func(terms.terms[i])
    end
    return sum
end
Base.:*(op::Operator, terms::Terms)=operator_multiplication_parallel_handler(op,terms,Val(length(terms)>30))

#Or the operator can be called directly
(op::Operator)(herm::Terms)::Terms=op*herm
(op::Operator)(herm::Hermite)::Terms=op*herm
(op::Operator)(herm)::Terms=op*herm

#You can also raise an operator to an integer power
Base.:^(op::Operator,expon::Int)=begin
    result=1
    for j in 1:expon
        result*=op
    end
    return result
end

##Algebra for operators
#For two operators
Base.:*(op1::Operator, op2::Operator)=Operator((herm::Hermite)->op1*(op2*herm))
Base.:+(op1::Operator, op2::Operator)=Operator((herm::Hermite)->(op1*herm+op2*herm))
Base.:-(op1::Operator, op2::Operator)=Operator((herm::Hermite)->op1*herm-(op2*herm))
Base.:/(op1::Operator, op2::Operator)=Operator((herm::Hermite)->op1*(1/(op2*herm)))

#For a non-operator and an operator
Base.:*(@nospecialize(factor),op::Operator)=Operator((herm::Hermite)->factor*(op*herm))
Base.:*(op::Operator,factor)=factor*op
Base.:+(@nospecialize(factor),op::Operator)=Operator((herm::Hermite)->factor+(op*herm))
Base.:+(op::Operator,factor)=factor+op
Base.:-(@nospecialize(factor),op::Operator)=Operator((herm::Hermite)->factor-(op*herm))
Base.:-(op::Operator,factor)=-(factor-op)
Base.:/(@nospecialize(factor),op::Operator)=Operator((herm::Hermite)->factor/(op*herm))
Base.:/(op::Operator,factor)=Operator((herm::Hermite)->(op*herm)/factor)

Base.:-(op::Operator)=Operator((herm::Hermite)->-(op*herm))

Symbolics.substitute(op::Operator, args...)=Operator((herm::Hermite)->Symbolics.substitute(op*herm, args...))

#Define the identity Operator
Base.one(::Operator)=Operator((herm::Hermite)->Terms(herm))
Base.one(::Type{Operator})=Operator((herm::Hermite)->Terms(herm))

#X raises the order of a Hermite basis state (like multiplying by x)
@variables x t
function Xfunc(herm::Hermite)::Terms
    Terms(Hermite(herm.factor, herm.a, herm.b, herm.order+1))
end
const X=Operator(Xfunc)
#Derivative with respet to x
function Dxfunc(herm::Hermite)::Terms
    herm.order==0 ? Terms([Hermite(herm.factor*herm.a*2,herm.a,herm.b,1)]) : Terms([Hermite(herm.factor*herm.order,herm.a,herm.b,herm.order-1),Hermite(herm.factor*herm.a*2,herm.a,herm.b,herm.order+1)])
end
const Dx=Operator(Dxfunc)
