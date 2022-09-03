# Generates a Maxima program for computing the integrals of the Lagrange bases,
# i.e., the weights for the Newton-Cotes formulas.

"""
Usage:
   julia> open("newton-cotes.mac", "w") do f newton_cotes(12, f) end
   (%i1) load("newton-cotes.mac");
   (%i2) bfloat(H);
"""
function newton_cotes(nmax, stream)
    for n in 1:nmax, i in 1:n+1
        r = i - 1
        print(stream, "H_$(n)_$i : (-1)**($n-$r) / ($(r)! * ($n-$r)!) * integrate(")
        for j in 0:n
            if r != j
                print(stream, "(t-$j)")
                if j != n && !(j == n - 1 && r == n)
                    print(stream, "*")
                end
            end
        end
        println(stream, ", t, 0, $n);")
    end
    for n in 1:nmax
        print(stream, "H_$n : [")
        for i in 1:n+1
            print(stream, "H_$(n)_$i")
            if i != n + 1
                print(stream, ", ")
            end
        end
        println(stream, "];")
    end
    print(stream, "H: [")
    for n in 1:nmax
        print(stream, "H_$n")
        if n != nmax
            print(stream, ", ")
        end
    end
    println(stream, "];")
    println(stream, "display2d: false;");
end
