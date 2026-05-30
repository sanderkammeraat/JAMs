

#Unify force input using abstract types for e.g.  pair forces ("K") and external forces ("M")
abstract type K end

struct variantK1<:K
    test
end

struct variantK2<:K
    test2
end

abstract type M end

struct variantM1<:M
end

combination = [variantK1(1), variantK2(1), variantM1()]

for variant in combination
    
    if isa(variant, K)
        println("$variant is of K")
    else
        println("$variant is not of K")
    end
end

Kvariants = [ Kvar for Kvar in combination if isa(Kvar, K)]