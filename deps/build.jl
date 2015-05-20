on_windows = @windows ? true : false

oldwdir = pwd()

pkgdir = Pkg.dir("Nemo") 

if on_windows
   wdir = "$pkgdir\\src"
   wdir2 = split(wdir, "\\")
   s = lowercase(shift!(wdir2)[1])
   unshift!(wdir2, string(s))
   unshift!(wdir2, "")
   wdir2 = join(wdir2, "/") 
else
   wdir = "$pkgdir/src"
end

cd(wdir)

if on_windows
   oldpth = ENV["PATH"]
   pth = split(oldpth, ";")
   shift!(pth)
   shift!(pth)
   pth = join(pth, ";")

   start = ENV["COMSPEC"]

   if Int == Int32
      ENV["MSYSTEM"]="MINGW32"
   else
      ENV["MSYSTEM"]="MINGW64"
   end
else
   # install M4
   run(`wget http://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.bz2`)
   run(`tar -xvf m4-1.4.17.tar.bz2`)
   run(`rm m4-1.4.17.tar.bz2`)
   cd("$wdir/m4-1.4.17")
   run(`./configure --prefix=$wdir`)
   run(`make`)
   run(`make install`)
   cd(wdir)
   run(`rm -rf m4-1.4.17`)
end

# install M4

if !on_windows
   run(`wget http://ftp.gnu.org/gnu/m4/m4-1.4.17.tar.bz2`)
   run(`tar -xvf m4-1.4.17.tar.bz2`)
   run(`rm m4-1.4.17.tar.bz2`)
   cd("$wdir/m4-1.4.17")
   run(`./configure --prefix=$wdir`)
   run(`make`)
   run(`make install`)
end

cd(wdir)

# install MPIR

run(`wget http://mpir.org/mpir-2.7.0-alpha11.tar.bz2`)
run(`tar -xvf mpir-2.7.0-alpha11.tar.bz2`)
run(`rm mpir-2.7.0-alpha11.tar.bz2`)
cd("$wdir/mpir-2.7.0")

if on_windows
   ENV["PATH"] = pth
   if Int == Int32
      run(`$start /e:4096 /c sh configure --prefix=$wdir2 --enable-gmpcompat --disable-static --enable-shared ABI=32`)
   else
      run(`$start /e:4096 /c sh configure --prefix=$wdir2 --enable-gmpcompat --disable-static --enable-shared ABI=64`)
   end
   run(`$start /e:4096 /c sh -c "make -j"`)
   run(`$start /e:4096 /c sh -c "make install"`) # naturally autotools fails to do this correctly
   run(`$start /e:4096 /c sh -c "cp .libs/libgmp*.dll $wdir2/lib"`)  # so we do it ourselves
   run(`$start /e:4096 /c sh -c "cp .libs/libmpir*.dll $wdir2/lib"`)  # so we do it ourselves
   ENV["PATH"] = oldpth
else
   run(`./configure --prefix=$wdir M4=$wdir/bin/m4 --enable-gmpcompat --disable-static --enable-shared`)
   run(`make -j4`)
   run(`make install`)
   cd(wdir)
   run(`rm -rf mpir-2.7.0`)
   run(`rm -rf bin`)
end

cd(wdir)

# install MPFR

run(`wget http://www.mpfr.org/mpfr-current/mpfr-3.1.2.tar.bz2`)
run(`tar -xvf mpfr-3.1.2.tar.bz2`)
run(`rm mpfr-3.1.2.tar.bz2`)
cd("$wdir/mpfr-3.1.2")

if on_windows
   ENV["PATH"] = pth
   run(`$start /e:4096 /c sh configure --prefix=$wdir2 --with-gmp=$wdir2 --disable-static --enable-shared`)
   run(`$start /e:4096 /c sh -c "make -j"`)
   run(`$start /e:4096 /c sh -c "make install"`) # naturally autotools fails to do this correctly
   run(`$start /e:4096 /c sh -c "cp src/.libs/libmpfr*.dll $wdir2/lib"`)  # so we do it ourselves
   ENV["PATH"] = oldpth
else
   run(`./configure --prefix=$wdir --with-gmp=$wdir --disable-static --enable-shared`)
   run(`make -j4`)
   run(`make install`)
   cd(wdir)
   run(`rm -rf mpfr-3.1.2`)
end

cd(wdir)

# install ANTIC

try
  run(`git clone https://github.com/wbhart/antic.git`)
except
  run(`cd antic ; git pull`)
end          

# install FLINT
try
  run(`git clone https://github.com/wbhart/flint2.git`)
except
  run(`cd flint2 ; git pull`)
end          

if on_windows
   cd("$wdir\\flint2")
   ENV["PATH"] = pth
   if Int == Int32
      run(`$start /e:4096 /c sh configure --extensions="$wdir2\antic" --prefix=$wdir2 --disable-static --enable-shared --with-mpir=$wdir2 --with-mpfr=$wdir2 ABI=32`)
   else
      run(`$start /e:4096 /c sh configure --extensions="$wdir2\antic" --prefix=$wdir2 --disable-static --enable-shared --with-mpir=$wdir2 --with-mpfr=$wdir2 ABI=64`)
   end
   run(`$start /e:4096 /c sh -c "make -j"`)
   run(`$start /e:4096 /c sh -c "make install"`)
   ENV["PATH"] = oldpth
else
   cd("$wdir/flint2")
   run(`./configure --prefix=$wdir --extensions="$wdir/antic" --disable-static --enable-shared --with-mpir=$wdir --with-mpfr=$wdir`)
   run(`make -j4`)
   run(`make install`)
end

cd(wdir)

# INSTALL ARB 

try
  run(`git clone -b julia https://github.com/thofma/arb.git`)
except
  run(`cd arb ; git pull`)
end          
 
if on_windows
else
   cd("$wdir/arb")
   run(`./configure --prefix=$wdir --disable-static --enable-shared --with-mpir=$wdir --with-mpfr=$wdir --with-flint=$wdir`)
   run(`make -j4`)
   run(`make install`)
   cd(wdir)
end

# install PARI

try
  run(`git clone http://pari.math.u-bordeaux.fr/git/pari.git`)
except
  run(`cd pari ; git pull`)
end  

if on_windows
else
   cd("$wdir/pari")
   env_copy = copy(ENV)
   env_copy["LD_LIBRARY_PATH"] = "$wdir/lib"
   config_str = `./Configure --prefix=$wdir --with-gmp=$wdir --mt=pthread`
   config_str = setenv(config_str, env_copy)
   run(config_str)
   run(`make -j4 gp`)
   run(`make doc`)
   run(`make install`)
end

cd(wdir)

if on_windows
   push!(Libdl.DL_LOAD_PATH, "$pkgdir\\src\\lib")
else
   push!(Libdl.DL_LOAD_PATH, "$pkgdir/src/lib")
end

cd(oldwdir)

