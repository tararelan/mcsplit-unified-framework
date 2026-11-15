import config as C
D=C.load("dal","lv").set_index(["a","b"]); M=C.load("mcsplit","lv").set_index(["a","b"])
pairs=[("g12","g45"),("g21","g66"),("g46","g47"),("g6","g45")]
for p in pairs:
    for a,b in [p]:
        key=[k for k in D.index if k[0].startswith(a) and k[1].startswith(b)]
        for k in key:
            print(k, "DAL uni", D.loc[k,"size"], "DAL ref", D.loc[k,"ref_size"],
                  "| McSplit uni", M.loc[k,"size"], "McSplit ref", M.loc[k,"ref_size"],
                  "| DAL nodes", D.loc[k,"nodes"], "McSplit nodes", M.loc[k,"nodes"])