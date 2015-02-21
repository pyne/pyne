from subprocess import Popen, PIPE, STDOUT

def test_alphad():
    args = ("./alphad", "alphad.inp")
    p = Popen(args, stdout=PIPE, stdin=PIPE, stderr=STDOUT)    
    grep_stdout = p.communicate(input=b'\n\n\n\n\n\n')[0]    
    # or this: grep_stdout = p.communicate(input=b'alphad.inp\n\n\n\n\n')[0]
    print(grep_stdout)

test_alphad()
