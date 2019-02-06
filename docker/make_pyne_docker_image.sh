#! /bin/bash

#_____________________________________________________________________________
# usage
# call if option -h, -help or --help, display usage and quit
function usage ()
{
echo "---------------------------------------------------------------------------"
echo "---------------------- make_pyne_docker_image usage: ----------------------"
echo " "
echo "  -h|-help|--help: print help/usage "
echo " "
echo "  --moab|-moab|moab : add MOAB to the build "
echo " "
echo "  --pymoab|-pymoab|pymoab : add pyMOAB to the build (also enable MOAB) "
echo " "
echo "  --dagmc|-dagmc|dagmc : add DAGMC to the build (also enable MOAB) "
echo " "
echo "  --deps-only|-deps-only|deps-only|--deps|-deps|deps : only build "\
    "dependencies (PyNE will not be build/installed)"
echo " "
echo "  -tag=*|--tag=*|tag=*|-t=*|--t=*|t=*: set the tag name of the "\
            "created docker image, if not provided tag name will be form from "\
            "the install option : 
            \"ubuntu_18.04(_pyne/_pyne-deps)(moab/pymoab/)(dagmc/) \""
echo " "
echo "---------------------------------------------------------------------------"
exit 418 
}


# Parsing arguments
MOAB=false
DAGMC=false
PYMOAB=false
PYNE=true

for arg in "$@"; do
    case $arg in
        -h|-help|--help )
            usage 
            ;;
        --moab|-moab|moab )
            MOAB=true
            ;;
        --dagmc|-dagmc|dagmc )
            DAGMC=true 
            MOAB=true 
            ;;
        --pymoab|-pymoab|pymoab )
            PYMOAB=true 
            MOAB=true 
            ;;
        --deps-only|-deps-only|deps-only )
            PYNE=false 
            ;;
        --deps|-deps|deps )
            PYNE=false 
            ;;
        -tag=*|--tag=*|tag=*|-t=*|--t=*|t=*)
            img_tag="${arg#*=}"
            shift # past argument=value
            ;;
    esac
done

# Parse option, form arguments
pyne_option='_pyne'
if  $MOAB; then
    moab_arg='--build-arg build_moab=YES'
    moab_option='_moab'
fi
if  $DAGMC; then
    dagmc_arg='--build-arg build_dagmc=YES'
    dagmc_option='_dagmc'
    moab_option=''
fi
if  $PYMOAB; then
    pymoab_arg='--build-arg enable_pymoab=YES'
    pymoab_option='_pymoab'
    moab_option-''
fi
if  ! $PYNE ; then
    pyne_arg='--build-arg build_pyne=NO'
    pyne_option='_pyne-deps'
fi

# if no image tag, create it
if [ -z "$img_tag" ]; then
    img_tag=ubuntu_18.04${pyne_option}${moab_option}${pymoab_option}${dagmc_option}
fi

docker build -t ${img_tag} -f ubuntu_18.04-dev.dockerfile ${moab_arg} ${dagmc_arg} ${pymoab_arg} ${pyne_arg} .



exit 0

EOF
##############################################################################

%DOC

# Exit status
In `bash`, a correct exectution of a script return an exit status equals to `0` (that why a *C*/*C++* program end by `return 0;`). There is NO standard for other exit status, so I (Josselin Massot) use the list of status code in HTTP for exit status.

* `404` : *Not Found*, a test find a file which doesn't exist.
* `418` : *I'm a teapot*, if only return help (I don't find a correct exit status for this).
* `501` : *Not Implemented*, test a feature which doesn't implemented on OS.
* `505` : *HTTP Version NOt supported*, the version of a library is NOt the one is expected.
* `507` : *Insufficient storage*, an error on making a directory, maybe because of the insufficient storage.

