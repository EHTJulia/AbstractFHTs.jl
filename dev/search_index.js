var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = AbstractFHTs","category":"page"},{"location":"#AbstractFHTs","page":"Home","title":"AbstractFHTs","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for AbstractFHTs.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [AbstractFHTs]","category":"page"},{"location":"#AbstractFHTs.fht","page":"Home","title":"AbstractFHTs.fht","text":"fht(A [, dims])\n\nPerforms a multidimensional Fast Hartley Transform (FHT) of the array A. The optional dims argument specifies an iterable subset of dimensions  (e.g. an integer, range, tuple, or array) to transform along. Most efficient if the size of A along the transformed dimensions is a product of small primes; see Base.nextprod. See also plan_fht() for even greater efficiency.\n\nA one-dimensional FHT computes the one-dimensional discrete Hartley transform (DHT) as defined by\n\noperatornameDFT(A)k =\n  sum_n=1^operatornamelength(A)\n  casleft(+ifrac2pi\n  (n-1)(k-1)operatornamelength(A) right) An\n\nwhere \\cas is the cosine-and-sine function, or alternatively called Hartley kernel, defined by\n\n\\cas(x) = \\cos(x) + \\sin(x).\n\nA multidimensional FHT simply performs this operation along each transformed dimension of A.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.fht!","page":"Home","title":"AbstractFHTs.fht!","text":"fht!(A [, dims])\n\nSame as fht, but operates in-place on A.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.ifht","page":"Home","title":"AbstractFHTs.ifht","text":"ifft(A [, dims])\n\nMultidimensional inverse FHT.\n\nA one-dimensional inverse FHT computes the one-dimensional inverse DHT (IDHT) defined by as\n\noperatornameIDHT(A)k = frac1operatornamelength(A) operatornameDHT(A)k\n\nA multidimensional inverse FHT simply performs this operation along each transformed dimension of A.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.ifht!","page":"Home","title":"AbstractFHTs.ifht!","text":"ifht!(A [, dims])\n\nSame as ifht, but operates in-place on A.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.plan_fht","page":"Home","title":"AbstractFHTs.plan_fht","text":"plan_fht(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)\n\nPre-plan an optimized FHT along given dimensions (dims) of arrays matching the shape and type of A.  (The first two arguments have the same meaning as for fht(@ref)).) Returns an object P which represents the linear operator computed by the FHT, and  which contains all of the information needed to compute fht(A, dims) quickly.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.plan_fht!","page":"Home","title":"AbstractFHTs.plan_fht!","text":"plan_fht!(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)\n\nSame as plan_fht, but operates in-place on A.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.plan_ifht","page":"Home","title":"AbstractFHTs.plan_ifht","text":"plan_ifht(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)\n\nSame as plan_fht, but produces a plan that performs inverse transforms ifht.\n\n\n\n\n\n","category":"function"},{"location":"#AbstractFHTs.plan_ifht!","page":"Home","title":"AbstractFHTs.plan_ifht!","text":"plan_ifht!(A [, dims]; flags=FFTW.ESTIMATE, timelimit=Inf)\n\nSame as plan_ifht, but operates in-place on A.\n\n\n\n\n\n","category":"function"},{"location":"#Base.size-Tuple{AbstractFHTs.DHTPlan, Any}","page":"Home","title":"Base.size","text":"size(p::DHTPlan, [dim])\n\nReturn the size of the input of a plan p, optionally at a specified dimenion dim.\n\n\n\n\n\n","category":"method"}]
}
