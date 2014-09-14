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
end

# install MPIR

run(`wget http://mpir.org/mpir-2.7.0-alpha10.tar.bz2`)
run(`tar -xvf mpir-2.7.0-alpha10.tar.bz2`)
run(`rm mpir-2.7.0-alpha10.tar.bz2`)
cd("$wdir/mpir-2.7.0")

if on_windows
   run(`sed -i 's/data.rel.ro.local,"aw",@progbits/data.rel.ro.local,"aw"/g' $wdir2/mpir-2.7.0/mpn/x86_64/x86_64-defs.m4`)
   ENV["PATH"] = pth
   if Int == Int32
      run(`$start /e:4096 /c sh configure --prefix=$wdir2 --enable-gmpcompat --disable-static --enable-shared ABI=32`)
   else
      run(`$start /e:4096 /c sh configure --prefix=$wdir2 --enable-gmpcompat --disable-static --enable-shared ABI=64`)
   end
   run(`$start /e:4096 /c sh -c "make -j"`)
   run(`$start /e:4096 /c sh -c "make install"`)
   ENV["PATH"] = oldpth
else
   run(`./configure --prefix=$wdir --enable-gmpcompat --disable-static --enable-shared`)
   run(`make -j4`)
   run(`make install`)
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
   run(`$start /e:4096 /c sh -c "make install"`)
   ENV["PATH"] = oldpth
else
   run(`./configure --prefix=$wdir --with-gmp=$wdir --disable-static --enable-shared`)
   run(`make -j4`)
   run(`make install`)
end

cd(wdir)

# install FLINT

run(`git clone https://github.com/wbhart/flint2.git`)

if on_windows
   cd("$wdir\\flint2")
   ENV["PATH"] = pth
   if Int == Int32
      run(`$start /e:4096 /c sh configure --prefix=$wdir2 --disable-static --enable-shared --with-mpir=$wdir2 --with-mpfr=$wdir2 ABI=32`)
   else
      run(`$start /e:4096 /c sh configure --prefix=$wdir2 --disable-static --enable-shared --with-mpir=$wdir2 --with-mpfr=$wdir2 ABI=64`)
   end
   run(`$start /e:4096 /c sh -c "make -j"`)
   run(`$start /e:4096 /c sh -c "make install"`)
   ENV["PATH"] = oldpth
else
   cd("$wdir/flint2")
   run(`./configure --prefix=$wdir --disable-static --enable-shared --with-mpir=$wdir --with-mpfr=$wdir`)
   run(`make -j4`)
   run(`make install`)
end

cd(wdir)

if on_windows
   push!(DL_LOAD_PATH, "$pkgdir\\src\\lib")
else
   push!(DL_LOAD_PATH, "$pkgdir/src/lib")
end

cd(oldwdir)

