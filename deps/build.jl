pkgdir = Pkg.dir("Nemo") 
pwd = "$pkgdir/src"

# install MPIR

run(`wget http://mpir.org/mpir-2.7.0-alpha10.tar.bz2`)
run(`tar -xvf mpir-2.7.0-alpha10.tar.bz2`)
run(`rm mpir-2.7.0-alpha10.tar.bz2`)
cd("$pwd/mpir-2.7.0")
run(`./configure --prefix=$pwd --enable-gmpcompat --disable-static --enable-shared`)
run(`make -j4`)
run(`make install`)
cd(pwd)

# install MPFR

run(`wget http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.bz2`)
run(`tar -xvf mpfr-3.1.2.tar.bz2`)
run(`rm mpfr-3.1.2.tar.bz2`)
cd("$pwd/mpfr-3.1.2")
run(`./configure --prefix=$pwd --with-gmp=$pwd --disable-static --enable-shared`)
run(`make -j4`)
run(`make install`)
cd(pwd)

# install FLINT

run(`git clone https://github.com/wbhart/flint2.git`)
cd("$pwd/flint2")
run(`./configure --prefix=$pwd --disable-static --enable-shared --with-mpir=$pwd --with-mpfr=$pwd`)
run(`make -j4`)
run(`make install`)
cd(pwd)

push!(DL_LOAD_PATH, "$pwd/lib")

