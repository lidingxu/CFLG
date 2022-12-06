from pathlib import Path
import os
import sys

# 1: solved. 2: unsolved
def main(argv):
    result_dir=argv[1] 
    instance=argv[2]
    algo=argv[3]
    cover=argv[4]
    result_path = result_dir + "/" + instance + "." + algo+ "." + cover
    #print(result_dir , instance , algo , cover)
    if os.path.isfile(result_path) and os.path.getsize(result_path) > 0:
        return 1
    else:
        if os.path.isfile(result_path):
            os.remove(result_path)
        return 2
    return 0

if __name__ == "__main__":
    rtcode = main(sys.argv)
    sys.exit(rtcode)