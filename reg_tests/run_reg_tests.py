# =============================================================================
# Standard Python modules                                           
# =============================================================================
import os, sys, argparse

# =============================================================================
# Extension modules
# =============================================================================
import mdo_regression_helper as reg

# define scripts to run:
test_scripts = ['test_curves.py','test_surfaces.py','test_volumes.py']
module_name = 'pyspline'

# Get the optional commandline arguments:
parser = argparse.ArgumentParser()
parser.add_argument("-mode",default='compare',choices=['train','compare'],
                    help='Train generates reference file. Compare runs test')
parser.add_argument("-diff_cmd",default='xxdiff',
                    help='Command to run for displaying diff. Default: xxdiff')
parser.add_argument("-nodiff", action='store_true', help='Suppress\
 displaying the comparison if not successful')
args = parser.parse_args()

mode = args.mode
diff_cmd = args.diff_cmd
nodiff = args.nodiff

if mode == 'train':
    try:
        os.remove('%s_reg.ref'%(module_name))
    except OSError:
        pass

    # Run each script
    for test in test_scripts:
        os.system('python %s >> %s_reg.ref'%(test, module_name))
    # end for
            
    # If we're training, we done (no comparison)
    sys.exit(0)
else:
    try:
        os.remove('%s_reg'%(module_name))
    except OSError:
        pass

    for test in test_scripts:
        os.system('python %s >> %s_reg 2>&1'%(test, module_name))
    # end for

    # Do the comparison (reference file must be first)
    res = reg.reg_file_comp('%s_reg.ref'%(module_name),'%s_reg'%(module_name))
# end if

# Set the proper return codes for the script running this:
if res == 0: #reg.REG_FILES_MATCH
    print('%s: Success!'%(module_name))
elif res == 1: #reg.REG_FILES_DO_NOT_MATCH
    print('%s: Failure!'%(module_name))
    if not nodiff:
        os.system('%s %s_reg.ref %s_reg'%(diff_cmd, module_name, module_name))
    else:
        os.system('cat %s_reg'%(module_name))
   
elif res == -1: #reg.REG_ERROR
    print('%s: Error in regression. Missing files.'%(module_name))
# end if

# Exit with code from reg_file_comp:
sys.exit(res)
