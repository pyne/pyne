from subprocess import Popen, PIPE, STDOUT

def test_delta():
    args = ("./delta")
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    grep_stdout = p.communicate(input=b'delta.dat\nout\n\n\n')[0]    
    # or this: grep_stdout = p.communicate(input=b'alphad.inp\n\n\n\n\n')[0]
    print(grep_stdout)

test_delta()
