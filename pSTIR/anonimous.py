def f(x, y):
    if len(x) > len(y):
        return 1
    else:
        return 0
def g(u):
    print(u('house'))

y = 'car'
g(lambda x: f(x,y))


