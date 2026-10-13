from pfac import fac

fac.SetAtom("Fe")

fac.Config(
    "1s2 2s2 2p6 3s2 3p6 3d6 4s1",
    group="FeII_ref"
)

fac.Config(
    "1s2 2s2 2p6 3s2 3p5 3d7 4s1",
    group="FeII_ref"
)

fac.Config(
    "1s2 2s2 2p6 3s2 3p5 3d6 4s2",
    group="FeII_ref"
)

fac.Config(
    "1s2 2s2 2p5 3s2 3p6 3d7 4s1",
    group="FeII_ref"
)

fac.Config(
    "1s2 2s2 2p5 3s2 3p6 3d6 4s2",
    group="FeII_ref"
)

fac.Config(
    "1s1 2s2 2p6 3s2 3p6 3d7 4s1",
    group="FeII_ref"
)

fac.Config(
    "1s1 2s2 2p6 3s2 3p6 3d6 4s2",
    group="FeII_ref"
)

fac.ConfigEnergy(0)

fac.OptimizeRadial(["FeII_ref"])

fac.ConfigEnergy(1)

fac.Structure("FeII_mu0.lev.b",["FeII_ref"])

fac.MemENTable("FeII_mu0.lev.b")
fac.PrintTable("FeII_mu0.lev.b","FeII_mu0.lev",1)

with open("FeII_mu0.lev", "r") as f:
    for line in f:
        if line.startswith("E0"):
            print("\nGround-state energy reference:")
            print(line.strip())
            break
