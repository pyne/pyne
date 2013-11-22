#build Stack Script
#need CUBIT Version#, optional pointer to cubit stuff
#CGM verion #, optional pointer to built stuff or tgz (or we can get from svn)
#hdf5 stuff
#MOAB pointer
#MCNP #
#place to put stuff

##################################
#
# CODE SPECIFIC BUILD INFO
#
##################################

# default variable values
no_mcnp=""
mcnp_patch=""
mcnp_source=$HOME/LANL/MCNP5/Source
mcnp_is_patched=""

code_long_options="no_mcnp,  
    mcnp_source:,
    mcnp_is_patched,
    mcnp_patch:"

code_help_msg="    --no_mcnp			Do not build MCNP
    --mcnp_source=SRC           MCNP Source Code located in SRC ["'$HOME'"/LANL/MCNP5/Source]      
    --mcnp_is_patched           MCNP Source already been patched
    --mcnp_patch=PATCH          Location of patch file PATCH
"
function parseCodeOption(){

    case $1 in
        --no_mcnp)		    no_mcnp="yes";;	
        --mcnp_is_patched) 	mcnp_is_patched="yes";;
        --mcnp_source)	    mcnp_source=$2 ; shift_arg=true;;
        --mcnp_patch)   	mcnp_patch=$2; shift_arg=true;;
        
        (*) code_option_found=false;
    esac
}

function buildCode() {
    if [ "$no_mcnp" != "" ]
    then
        echo "Quitting without building mcnp"
	    exit
    fi
    
    if [ -z $mcnp_patch ]
    then
        mcnp_patch=$dagmc_prefix/MCNP5/patch/dagmc.patch.5.1.60
        echoErr "Set patch to "$mcnp_patch
    fi
    
#patch mcnp5
    
    cd $mcnp_source
    if [ -z $mcnp_is_patched ]
    then
        echoErr "applying MCNP5 patch"
	    patch -p 1 < $mcnp_patch
    fi
    
    make build CONFIG="seq plot gfortran dagmc" FC=gfortran MARCH=M64 \
        MOAB_DIR=$moab_prefix CUBIT_DIR=$cubit_bin_path \
        DAGMC_DIR=$dagmc_prefix/MCNP5/dagmc
    
    mkdir -p $prefix/bin
    cp $mcnp_source/../bin/mcnp5 $prefix/bin/dag-mcnp5
}

###############################
#
# COMMON FOR ALL CODES
#
###############################

#default variable values
prefix=$HOME/dagmc_bld
clean="no"
make_parallel=""
_V=0

path_to_cubit=""

enable_static="no"

cgm_repo=https://svn.mcs.anl.gov/repos/ITAPS/cgm/

cgm_code="none"
cgm_special="none"
cgm_revision=""
cgm_branch=""
cgm_tag=""


moab_repo=https://bitbucket.org/fathomteam/moab.git
moab_code="none"
moab_special="none"
moab_revision=""
moab_branch=""
moab_tag=""

dagmc_code="none"
dagmc_special="none"
dagmc_tag=""

echoErr() { echo "$@" 1>&2; }

checkErr() {
    if [ $? -gt 0 ]; then
        echoErr "ERROR! $@"
        exit    
    fi
}

# function to unpack a compressed tarball and create a simlink to the src
function unTar(){   #param is file name
		cp $1 .
        checkErr "copying tgz failed"
        tgz=$(echo $1 | awk '{n = split($0,a,"/");  print a[n];}');
        echoErr "unzipping $tgz"
		tar -xzf $tgz
        checkErr "untarring tgz failed"
		tgz_dir=$(echo $tgz | awk '{split($0,a,".tar.gz");  print a[1];}');
		tgz_dir=$(echo $tgz_dir | awk '{split($0,a,".tgz");  print a[1];}');
        echoErr $tgz_dir
        ln -s $tgz_dir src
}

function svnCO(){   #param 1 is repo #param2 is code $3 is special
     if [ $2 == "revision" ]
        then              
            svn_dir=trunk 
            echoErr "Checking out rev $3 from $1/$svn_dir"
            svn co -r $3 $1/trunk
        elif [ $2 == "branch" ]
        then
            svn_dir=$3
            echoErr "Checking out from $1/branches/$svn_dir"
            svn co $1/branches/$svn_dir
        elif [ $2 == "tag" ]
        then 
            svn_dir=$3 
            echoErr "Checking out from $1/tags/$svn_dir"
            svn co $1/tags/$svn_dir           
        else
            svn_dir=trunk 
            echoErr "Checking out from $1/$svn_dir"
            svn co $1/$svn_dir
        fi	

        cd $svn_dir
		autoreconf -fi
		cd ..
		ln -fs $svn_dir src
}


function gitCO(){   #param 1 is repo url #param2 is code $3 is branch/name,4 is repo name
        git_cmd=""
        if [ $2 == "branch" ] 
        then
            git_cmd="-b $3"
        fi
        git clone $git_cmd $1
        cd $4
	if [ $2 == "tag" ] 
	then
	    git checkout tags/$3
	fi
        autoreconf -fi
        cd ..
        ln -fs $4 src
}

function buildCGM(){
    
    if [ $enable_static == "no" ]
    then
        if [ -z $cgm_installed ] && [ "$cgm_installed_dir" == "" ]
        then 
            
            
	        mkdir -p $cgm_prefix/bld
	        cd $cgm_prefix
            
	        if [ -z $cgm_tgz ]
	        then
                svnCO $cgm_repo $cgm_code $cgm_special
	        else
                unTar $cgm_tgz
	        fi
            
            checkErr "CGM setup failed"
            
	        cd bld		
	        ../src/configure --enable-optimize --disable-debug \
	            --with-cubit=$path_to_cubit/  \
	            --prefix=$cgm_prefix  
            
	        checkErr "CGM configure failed"
            make $make_parallel
            checkErr "CGM make failed"
            
	        make install 
            checkErr "CGM local install failed"
            
            echoErr "Successfully installed cgm"
	        cd $prefix
            
        else
	        if [ "$cgm_installed_dir" != "" ]
	        then
                echoErr $cgm_installed_dir
		        cgm_prefix=$cgm_installed_dir
	        fi
        fi
    fi
    
}


function buildHDF5() {
    if [ -z $hdf5_installed ] && [ "$hdf5_installed_dir" == "" ]
    then
	    mkdir -p $hdf5_prefix/bld
	    cd $hdf5_prefix
        
        unTar $hdf5_tgz
        
	    cd bld
	    ../src/configure --prefix=$hdf5_prefix	
        
        checkErr "HDF5 configure failed"
        
	    make $make_parallel
        checkErr "HDF5 make failed"
	    make install
        checkErr "HDF5 install failed"
	    cd $prefix
    else
	    if [ "$hdf5_installed_dir" != "" ]
	    then
		    hdf5_prefix=$hdf5_installed_dir
	    fi
    fi
}

function buildMOAB() {
    if [ -z $moab_installed ] && [ "$moab_installed_dir" == "" ]
    then 
        
        mkdir -p $moab_prefix/bld
	    cd $moab_prefix
        
	    if [ -z $moab_tgz ]	
        then 
            gitCO $moab_repo $moab_code $moab_special moab
	    else
	        unTar $moab_tgz
	    fi
	    cd bld
        
        if [ $enable_static == "no" ]
        then
       	    ../src/configure --enable-optimize --disable-debug \
	            --with-cgm=$cgm_prefix  \
	            --prefix=$moab_prefix     
        else
    	    ../src/configure --enable-optimize --disable-debug \
	            --enable-static  \
	            --prefix=$moab_prefix
        fi
        
        checkErr "MOAB Configure failed"
        
	    make $make_parallel
        checkErr "MOAB make failed"
	    make install
        checkErr "MOAB install failed"
        
        echoErr "MOAB Successfully Installed"
	    cd $prefix
    else
        
	    if [ "$moab_installed_dir" != "" ]
	    then
		    moab_prefix=$moab_installed_dir
	    fi
    fi
}

function buildDAGMC(){

#get DAGMC
    dagmc_repo=https://github.com/svalinn/DAGMC.git
    if [ -z $dagmc_installed ] && [ "$dagmc_installed_dir" == "" ]
    then 
        
	    cd $prefix
        
	    if [ -z $dagmc_tgz ]	
        then 
            gitCO $dagmc_repo $dagmc_code $dagmc_special DAGMC
	    else
	        unTar $dagmc_tgz
	    fi
        
        echoErr "DAGMC Successfully Installed"
    else
        
	    if [ "$dagmc_installed_dir" != "" ]
	    then
		    dagmc_prefix=$dagmc_installed_dir
	    fi
    fi
}

# options may be followed by one colon to indicate they have a required argument
if ! options=$(getopt -u -o  abhp:j: -l "
    help,
    prefix:,
    verbose,   
    enable_static,
    clean,

    path_to_cubit:,

    cgm_tgz:,
    cgm_installed,
    cgm_installed_dir:,
    cgm_revision:,
    cgm_branch:,    
    cgm_tag:,

    hdf5_tgz:,
    hdf5_installed,
    hdf5_installed_dir:,

    moab_tgz:,
    moab_installed,
    moab_installed_dir:,
    moab_branch:,
    moab_tag:,

    dagmc_tgz:,
    dagmc_installed,
    dagmc_installed_dir:,
    dagmc_branch:,
    dagmc_tag:

    $code_long_options"  -- "$@")
then
    # something went wrong, getopt will put out an error message
    exit 
fi

set -- $options

while [ $# -gt 0 ]
do
    code_option_found=true
    shift_arg=false

    parseCodeOption $1 $2

    if [ $code_option_found == false ]
    then
        case $1 in
            -h|--help)             help_flag="yes" ;;
            --prefix)              prefix=$2 ; shift;;
            --verbose)             _V=1;;
            --enable_static)       enable_static="yes";;
            --clean)               clean="yes";;
            -j)                    make_parallel="-j $2"; shift;;
            
            --cgm_tgz)		       cgm_tgz=$2; shift;;
            --cgm_installed)       cgm_installed="yes";;
            --cgm_installed_dir)   cgm_installed_dir=$2;shift;;
            --cgm_revision)		   cgm_code="revision";cgm_special=$2; shift;;
            --cgm_branch)          cgm_code="branch";cgm_special=$2; shift;;
            --cgm_tag)             cgm_code="tag";cgm_special=$2; shift;;
            
            --path_to_cubit)       path_to_cubit=$2;cubit_bin_path=$path_to_cubit/bin;shift;;
            
            --hdf5_tgz)		       hdf5_tgz=$2; shift;;
            --hdf5_installed)      hdf5_installed="yes";;
            --hdf5_installed_dir)  hdf5_installed_dir=$2;shift;;
            
            --moab_tgz)		       moab_tgz=$2; shift;;
            --moab_installed)      moab_installed="yes";;
            --moab_installed_dir)  moab_installed_dir=$2;shift;;
            --moab_branch)         moab_code="branch";moab_special=$2; shift;;
            --moab_tag)            moab_code="tag";moab_special=$2; shift;;
            
            
            --dagmc_tgz)		   dagmc_tgz=$2; shift;;
            --dagmc_installed)     dagmc_installed="yes";;
            --dagmc_installed_dir) dagmc_installed_dir=$2;shift;;
            --dagmc_branch)        dagmc_code="branch";dagmc_special=$2; shift;;
            --dagmc_tag)           dagmc_code="tag";dagmc_special=$2; shift;;
            
            (--) shift; break;;
            (-*) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
            (*) break;;
        esac
    fi

    if [ $shift_arg == true ]
    then
        shift
    fi

    shift

done

cubit_bin_path=$path_to_cubit/bin
cgm_prefix=$prefix/cgm/
hdf5_prefix=$prefix/HDF5
moab_prefix=$prefix/MOAB/
dagmc_prefix=$prefix/DAGMC

if [ "$help_flag" == "yes" ]
then
    echo "
Usage:    build_dagmc_stack.bash [OPTIONS]

        automatic building of DAGMCNP stack.
        Builds cgm, hdf5, MOAB and patched MCNP5 from svn download (default for
        cgm, MOAB) or from tgz (default hdf5)

Defaults for options are specified in brackets

cgm  repo : https://svn.mcs.anl.gov/repos/ITAPS/cgm/
hdf5 repo : None
moab repo : https://bitbucket.org/fathomteam/moab.git

OPTIONS:
    -h, --help                  display this help and exit 
    --prefix=PREFIX             install files in directory PREFIX ["'$HOME'"/dagmc_bld]
    --verbose                   get all svn and make info
    --clean                     Delete PREFIX before installing anything
    -j NUM, --j NUM             Use NUM processes for make
    --enable_static             Static build of MOAB [no] 

    --path_to_cubit             Path to Cubit

    --cgm_tgz=CGM               path to tgz file CGM		
    --cgm_installed             CGM already installed
    --cgm_installed_dir=CGM     cgm already installed in dir CGM
    --cgm_revision=REV          revision to pull from repo
    --cgm_branch=BRANCH         branch to pull from repo
    --cgm_tag=TAG               tag to pull from repo

    --hdf5_tgz=HDF5             path to tgz file HDF5		
    --hdf5_installed            HDF5 already installed
    --hdf5_installed_dir=HDF5   HDF5 already installed in dir HDF5  

    --moab_tgz=MOAB             path to tgz file MOAB
    --moab_installed            MOAB already installed
    --moab_installed_dir=MOAB   MOAB already installed in dir MOAB   
    --moab_branch=BRANCH        branch to pull from repo
    --moab_tag=TAG              tag to pull from repo

    --dagmc_tgz=DAGMC           path to tgz file DAGMC
    --dagmc_installed           MOAB already installed
    --dagmc_installed_dir=DAGMC MOAB already installed in dir DAGMC
    --dagmc_branch=BRANCH       branch to pull from repo
    --dagmc_tag=TAG             tag to pull from repo

$code_help_msg"     
    exit
fi

#set up verbose
if (($_V)); then
  echo "verbose"
  set -x
else
  exec 1>/dev/null
fi

if [ $clean == "yes" ]
then
    echoErr "Cleaning $prefix"
    rm -rf $prefix
    echoErr "$prefix cleaned"
fi

mkdir -p $prefix
cd $prefix

buildCGM; buildHDF5; buildMOAB; buildDAGMC

buildCode
