import re
def comment_remover(text):
    def replacer(match):
        s = match.group(0)
        if s.startswith('/'):
            return " " # note: a space and not an empty string
        else:
            return s
    pattern = re.compile(
        r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',
        re.DOTALL | re.MULTILINE
    )
    return re.sub(pattern, replacer, text)


f = open("define.h", "r")
txt = comment_remover(f.read())
text = "\""+"-".join([ll.rstrip() for ll in txt.splitlines() if ll.strip()])+"\""
f = open("define.dat","w")
f.write(text)
