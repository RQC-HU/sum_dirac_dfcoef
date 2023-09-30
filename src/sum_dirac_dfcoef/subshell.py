class Subshell:
    subshell_order = "spdfghiklmnoqrtuvwxyz"
    gto_label_order: "list[list[str]]" = []

    def lmnval(self, idx: int, nfun: int, istep: "list[int]", mval: "list[int]", nval: "list[int]") -> "tuple[list[int], list[int], list[int]]":
        # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/abacus/hergam.F#L218-236
        ix: "list[int]" = []
        iy: "list[int]" = []
        iz: "list[int]" = []
        for i in range(nfun):
            ix.append(idx - istep[i])
            iy.append(mval[i])
            iz.append(nval[i])
        return ix, iy, iz

    def carpow(self) -> "tuple[list[int], list[int], list[int]]":
        # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/abacus/hergam.F#L162-216
        istep, mval, nval = [], [], []
        for i in range(1, len(self.subshell_order) + 1):
            for j in range(1, i + 1):
                istep.append(i)
                mval.append(i - j)
                nval.append(j - 1)
        return istep, mval, nval

    def __init__(self):
        istep, mval, nval = self.carpow()
        # https://gitlab.com/dirac/dirac/-/blob/b10f505a6f00c29a062f5cad70ca156e72e012d7/src/abacus/herrdn.F#L4995-5021
        for idx, subshell in enumerate(self.subshell_order):
            if subshell == "s":
                self.gto_label_order.append(["s"])
            elif subshell == "p":
                self.gto_label_order.append(["px", "py", "pz"])
            elif subshell == "d":
                self.gto_label_order.append(["dxy", "dxy", "dxz", "dyy", "dyz", "dzz"])
            elif subshell == "f":
                self.gto_label_order.append(["fxxx", "fxxy", "fxxz", "fxyy", "fxyz", "fxzz", "fyyy", "fyyz", "fyzz", "fzzz"])
            else:
                nfun = (idx + 1) * (idx + 2) // 2
                ix, iy, iz = self.lmnval(idx + 1, nfun, istep, mval, nval)
                li = []
                for k in range(nfun):
                    li.append(self.subshell_order[idx] + str(ix[k]) + str(iy[k]) + str(iz[k]))
                self.gto_label_order.append(li)


subshell_order = Subshell()
