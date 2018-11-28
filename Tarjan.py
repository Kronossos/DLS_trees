def MakeSet(x):
    x.parent = x
    x.rank = 1


def TarjanOLCA(u):
    MakeSet(u)
    u.ancestor = u
    for v in u:

        if v is not u:
            TarjanOLCA(v)
        Union(u, v)
        Find(u).ancestor = u

    u.color = "black"
    # for v in P:
    #     if v.color == "black":
    #         print("Tarjan's Lowest Common Ancestor of " + u + " and " + v + " is " + Find(v).ancestor + ".")


def Find(x):
    if x.parent != x:
        x.parent = Find(x.parent)
    return x.parent


def Union(x, y):
    xRoot = Find(x)
    yRoot = Find(y)
    if xRoot.rank > yRoot.rank:
        yRoot.parent = xRoot
    elif xRoot.rank < yRoot.rank:
        xRoot.parent = yRoot
    else:
        yRoot.parent = xRoot
        xRoot.rank = xRoot.rank + 1
