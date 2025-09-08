using ASEconvert

ase.io.write("C.extxyz", ase.build.bulk("C", "diamond"))
ase.io.write("Si.extxyz", ase.build.bulk("Si", "diamond"))
ase.io.write("CsCl.extxyz", ase.build.bulk("CsCl", "cesiumchloride"; a=4.14))
