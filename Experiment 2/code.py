def add(a, b):
    return a + b

def sub(a, b):
    return a - b

def mul(a, b):
    return a * b 

def div(a, b):
    return a/b

x = int(input("Enter First number: "))
operation = input("Enter the operation: ")
y = int(input("Enter The second number "))

if operation == "+":
    ans = add(x, y)
    print(f"The addition is {ans}")
