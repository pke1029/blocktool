import numpy as np
import matplotlib.pyplot as plt

class Vert:

    def __init__(self, x, y, z=0, vid=None):
        self.x = x
        self.y = y
        self.z = z
        self.id = vid

    def set_id(self, vid):
        self.id = vid

    def __str__(self):
        return f"({self.x} {self.y} {self.z})\n"
        

class Edge:

    def __init__(self, v0, v1, x, y, z=0):
        self.v0 = v0
        self.v1 = v1
        self.x = x
        self.y = y
        self.z = z

    def __str__(self):
        return f"arc {self.v0.id} {self.v1.id} ({self.x} {self.y} {self.z})\n"

class Block:
    
    def __init__(self, v0, v1, v2, v3, v4, v5, v6, v7, nx, ny, nz, gx=1, gy=1, gz=1):
        self.verts = (v0, v1, v2, v3, v4, v5, v6, v7)
        # number of cells in each direction
        self.nx = nx
        self.ny = ny
        self.nz = nz 
        # cell expansion ratio
        self.gx = gx
        self.gy = gy
        self.gz = gz 

    def draw(self):
        x = [self.verts[0].x, self.verts[1].x, self.verts[2].x, self.verts[3].x, self.verts[0].x]
        y = [self.verts[0].y, self.verts[1].y, self.verts[2].y, self.verts[3].y, self.verts[0].y]
        plt.plot(x, y)
        plt.annotate(str(self.nx), ((x[0]+x[1])/2, (y[0]+y[1])/2), ha="center", va="center")
        plt.annotate(str(self.ny), ((x[1]+x[2])/2, (y[1]+y[2])/2), ha="center", va="center")

    def __str__(self):
        vids = " ".join(str(v.id) for v in self.verts)
        return f"hex ({vids}) ({self.nx} {self.ny} {self.nz}) simpleGrading ({self.gx} {self.gy} {self.gz})\n" 

class Face:

    def __init__(self, v0, v1, v2, v3):
        self.v0 = v0 
        self.v1 = v1 
        self.v2 = v2 
        self.v3 = v3

    def __str__(self):
        return f"({self.v0.id} {self.v1.id} {self.v2.id} {self.v3.id})\n"

class Patch:

    def __init__(self, name, ptype="patch", faces=[]):
        self.name = name
        self.ptype = ptype 
        self.faces = faces or []

    def add_face(self, face):
        self.faces.append(face)

    def __str__(self):
        text = f"{self.name}\n{{\n"
        text += f"\ttype {self.ptype};\n"
        text += f"\tfaces\n\t(\n"
        for f in self.faces:
            text += f"\t\t{f}"
        text += "\t);\n}\n\n"
        return text

class Patchpair:

    def __init__(self, patch1, patch2):
        self.patch1 = patch1 
        self.patch2 = patch2
    
    def __str__(self):
        return f"({self.patch1.name} {self.patch2.name})"

class Mesh:

    def __init__(self, zmin=0, zmax=1, nz=1, gz=1):
        self.zmin = zmin
        self.zmax = zmax
        self.nz = nz 
        self.gz = gz
        self.verts = []
        self.edges = []
        self.blocks = []
        self.patches = {}
        self.add_patch("top")
        self.add_patch("bottom", "wall")
        self.patchpairs = []

    def add_vert(self, x, y):
        v = Vert(x, y, self.zmin)
        w = Vert(x, y, self.zmax)
        n = len(self.verts)
        v.set_id(n)
        w.set_id(n+1)
        self.verts.append(v)
        self.verts.append(w)
        return v

    def set_vert_pos(self, vert, x, y):
        vid = vert.id 
        self.verts[vid].x = x
        self.verts[vid].y = y
        self.verts[vid+1].x = x
        self.verts[vid+1].y = y

    def add_edge(self, v0, v1, x, y):
        w0 = self.verts[v0.id + 1]
        w1 = self.verts[v1.id + 1]
        self.edges.append(Edge(v0, v1, x, y, self.zmin))
        self.edges.append(Edge(w0, w1, x, y, self.zmax))

    def add_block(self, v0, v1, v2, v3, nx=4, ny=4, gx=1, gy=1, dx=None, dy=None):
        if dx != None:
            nx = int(max(np.abs(v1.x-v0.x), np.abs(v3.x-v2.x)) // (dx-0.000001))
        if dy != None:
            ny = int(max(np.abs(v3.y-v0.y), np.abs(v2.y-v1.y)) // (dy-0.000001))
        w0 = self.verts[v0.id + 1]
        w1 = self.verts[v1.id + 1]
        w2 = self.verts[v2.id + 1]
        w3 = self.verts[v3.id + 1]
        self.blocks.append(Block(v0, v1, v2, v3, w0, w1, w2, w3, nx, ny, self.nz, gx, gy, self.gz))
        self.patches["bottom"].add_face(Face(v0, v1, v2, v3))
        self.patches["top"].add_face(Face(w0, w1, w2, w3))

    def add_patch(self, name, ptype="patch", verts=[]):
        faces = []
        for i in range(len(verts)-1):
            v0 = verts[i]
            v1 = verts[i+1]
            w0 = self.verts[v0.id+1]
            w1 = self.verts[v1.id+1]
            faces.append(Face(v0, v1, w1, w0))
        if name in self.patches.keys():
            self.patches[name].faces.extend(faces)
        else:
            self.patches[name] = Patch(name, ptype, faces)
    
    def add_patchpair(self, patch1, patch2):
        pp = Patchpair(patch1, patch2)
        self.patchpairs.append(pp)

    def new_box(self, xmin, ymin, xmax, ymax, pname0, pname1, pname2, pname3, dx=0.02, dy=0.02, gx=1, gy=1):
        v0 = self.add_vert(xmin, ymin)
        v1 = self.add_vert(xmax, ymin)
        v2 = self.add_vert(xmax, ymax)
        v3 = self.add_vert(xmin, ymax)
        nx = int((xmax-xmin) // (dx-0.000001))
        ny = int((ymax-ymin) // (dy-0.000001))
        self.add_block(v0, v1, v2, v3, nx, ny, gx, gy)
        self.add_patch(pname0, verts=[v0, v1])
        self.add_patch(pname1, verts=[v1, v2])
        self.add_patch(pname2, verts=[v2, v3])
        self.add_patch(pname3, verts=[v3, v0])

    def get_cell_count(self):
        return sum(b.nx * b.ny + b.nz for b in self.blocks)
    
    def write_section(self, title, items):
        body = "".join(f"\t{i}\n" for i in items)
        return f"{title}\n(\n{body});\n\n"

    def write_verts(self):
        return self.write_section("vertices", self.verts)

    def write_edges(self):
        return self.write_section("edges", self.edges)

    def write_blocks(self):
        return self.write_section("blocks", self.blocks)

    def write_patches(self):
        text = "boundary\n(\n"
        for p in self.patches.values():
            text += f"{p}"
        text += ");\n\n"
        return text
    
    def write_patchpairs(self):
        return self.write_section("mergePatchPairs", self.patchpairs)

    def write_all(self):
        header = """
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2412                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version         2.0;
    format          ascii;
    class           dictionary;
    object          blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale   1;

"""
        text = header
        text += self.write_verts()
        text += self.write_edges()
        text += self.write_blocks()
        text += self.write_patches()
        text += self.write_patchpairs()
        return text
    
    def show_verts(self):
        for v in self.verts:
            if v.z == self.zmin:
                plt.plot(v.x, v.y, ".r")
            elif v.z == self.zmax:
                plt.plot(v.x, v.y, ".b")
        plt.axis("equal")
        plt.show()
    
    def show_blocks(self):
        for b in self.blocks:
            b.draw()
        plt.axis("equal")
        plt.show()
