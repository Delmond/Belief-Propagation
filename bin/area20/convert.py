
f = open('factors.txt')

text = f.read()
lines = text.split('\n')

for i in range(len(lines)):
    lines[i] = lines[i].split()
    print(lines[i])

lines = lines[:-1]
for i in range(len(lines)):
    lines[i][0], lines[i][1] = lines[i][1], lines[i][0]

for i in range(len(lines)):
    lines[i] = '\t'.join(lines[i])

text = '\n'.join(lines)

print(text)
f2 = open('factors_csv.txt', 'w')
f2.write(text)