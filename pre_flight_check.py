import subprocess
import argparse
import os

parser = argparse.ArgumentParser(
                    prog = 'pre_flight_check',
                    description = 'check that the CircuitSeq tools are available',
                    epilog = 'please add a github issue for anything you run into!')

parser.add_argument('-s', '--sample_sheet')
args = parser.parse_args()

# https://stackoverflow.com/questions/2330245/python-change-text-color-in-shell
def hilite(string, status, bold):
    attr = []
    if status:
        # green
        attr.append('32')
    else:
        # red
        attr.append('31')
    if bold:
        attr.append('1')
    return '\x1b[%sm%s\x1b[0m' % (';'.join(attr), string)

# check the GPU status with NVIDIA's SMI command 
command = "nvidia-smi"

# Get the return status from the completed process
return_code = subprocess.call(command, shell=True)

# Check the return status
if return_code == 0:
    print("nvidia-smi ran successfully")
else:
    print("nvidia-smi failed with return code: " + str(return_code))
    raise NameError("Invalid return code for NVIDIA's SMI command (nvidia-smi)")

# now check that the sample sheet conforms to our standard
print(hilite("Checking sample sheet...",True,True))

def check_string(string):
    if string.isalnum() or "_" in string:
        return True
    else:
        return False

def file_exists(filename):
    return(os.path.isfile(filename))


sample_sheet = open(args.sample_sheet)
header = sample_sheet.readline().strip()
if (header != "position\tsample\treference"):
    print(hilite('We expect three column names, we saw: ' + str(header) + ', we expect: position<tab>sample<tab>reference<endline>',False,True))
    raise NameError("Invalid sample sheet")


for index,line in enumerate(sample_sheet):
    line_split = line.strip().split("\t")
    line = line.strip()
    if len(line_split) != 3:
        
        print(hilite('We expect three columns, we saw: ' + str(len(line_split)) + ', columns for line: ' + line + ' (line number ' + str(index + 2) + ')',False,True))
        raise NameError("Invalid sample sheet")
    
    if len(line_split[0]) != 2:
        print(hilite("We expect two-digit (i.e. 98, 09, 32) numbers in the first column, if less than 10 it should be padded with a zero, we saw: " + str(line_split[0]) + ", on line: " + line + " (line number " + str(index + 2) + ")",False,True))
        raise NameError("Invalid sample sheet")
    
    try:
        int_val = int(line_split[0])
    except:
        print(hilite("Unable to convert column 'position' into a number for line: " + line + " (line number " + str(index + 2) + ")",False,True))
        raise NameError("Invalid sample sheet")
    
    if not check_string(line_split[1]):
        print(hilite("We expect only letters, numbers, and underscores for samples, we saw: " + str(line_split[1]) + ", on line: " + line + " (line number " + str(index + 2) + ")",False,True))
        raise NameError("Invalid sample sheet")
    
    if not file_exists(line_split[2]):
        print(hilite("Your reference doesn't exist: " + str(line_split[2]) + ", on line: " + line + " (line number " + str(index + 2) + ")",False,True))
        raise NameError("Invalid sample sheet")
    

    
print(hilite("Done checking sample sheet! Success!",True,True))
