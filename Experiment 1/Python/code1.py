def add(a, b):
    return a + b
def sub(a, b):
    return a - b
def mul(a, b):
    return a * b
def div(a, b):
    return a / b

x = int(input(("Enter the first number: ")))
ope = input("Enter the operation: ")
y = int(input(("Enter the second number: ")))

if ope == '+':
    ans = add(x, y)
    print(ans)
elif ope == '-':
    ans = sub(x, y)
    print(ans)
elif ope == '*':
    ans = mul(x, y)
    print(ans)
elif ope == '/':
    ans = div(x, y)
    print("The answer is: ", ans)
else:
    print("Wrong operation entered")