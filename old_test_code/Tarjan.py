def MakeSet(x):
    x.parent = x
    x.rank = 1


def reference_search(u, P):
    for v in P:
        if u not in v: continue
        x = (v.index(u) + 1) % 2

        v_new = v[x]
        if v_new.color == "black":
            print("Tarjan's Lowest Common Ancestor of " + str(u.label()) + " and " + str(v_new.label()) + " is " + str(
                Find(v_new).ancestor.label()) + ".")


def label_search(u, P, label_dict):
    print("######", u.label(), "######")
    for v in P:

        left = v[0]
        right = v[1]

        # print(u.label())
        print("LEWY: ", left.label())
        print("PRWY: ", right.label())

        if len(set(u.label()) - set(left.label())) == 0:
            print("mam lewy")
            v_new = v[1]
        elif len(set(u.label()) - set(right.label())) == 0:
            print("mam prawy")
            v_new = v[0]
        else:
            continue



        if label_dict[v_new.label()].color == "black":
            print("Tarjan's Lowest Common Ancestor of " + str(u.label()) + " and " + str(v_new.label()) + " is " + str(
                Find(v_new).ancestor.label()) + ".")
        else:
            print("!")
            print()
            print()


def TarjanOLCA(u, P, label_dict, search_type=label_search):
    MakeSet(u)
    u.ancestor = u

    if not u.is_leaf():
        for v in (u.son("L"), u.son("R")):
            if v is not u:
                TarjanOLCA(v, P,label_dict)
                Union(u, v)
                Find(u).ancestor = u

    u.color = "black"

    search_type(u, P, label_dict)


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
