#%Module1.0#####################################################################
##
## MRO-3.5 modulefile
##
##
proc ModulesHelp { } {
        global dotversion

        puts stderr "\tAdds '/home/dongjiaoyang/miniconda3/envs/OESingleCell/bin' to your PATH environment variable"
        puts stderr "\tto your PATH environment variable.  This allows you to"
        puts stderr "\trun executables from the /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin."
}

module-whatis   "adds '/home/dongjiaoyang/miniconda3/envs/OESingleCell/bin' to your PATH environment variable"

prepend-path    PATH    /home/dongjiaoyang/miniconda3/envs/OESingleCell/bin

prepend-path -d " " CFLAGS "-I/home/dongjiaoyang/miniconda3/envs/OESingleCell/include"
prepend-path -d " " CPPFLAGS "-I/home/dongjiaoyang/miniconda3/envs/OESingleCell/include"

prepend-path -d " " LDFLAGS "-L/home/dongjiaoyang/miniconda3/envs/OESingleCell/lib"
prepend-path    LD_RUN_PATH /home/dongjiaoyang/miniconda3/envs/OESingleCell/lib
prepend-path    LD_LIBRARY_PATH /home/dongjiaoyang/miniconda3/envs/OESingleCell/lib
prepend-path    PYTHONPATH /home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/python3.7/site-packages
prepend-path    PATH /home/dongjiaoyang/miniconda3/envs/OESingleCell/lib/python3.7/site-packages/bin
