using NeuroAnalysis,Statistics,DataFrames,DataFramesMeta,Mmap,Images,StatsBase,Interact,CSV,MAT,DataStructures

using PyPlot, PyCall
@pyimport matplotlib.patches as patch
@pyimport matplotlib.pyplot as pylt
@pyimport descartes as desca
@pyimport matplotlib.colors as plcolor
@pyimport shapely as spl


LinearSegmentedColormap, ListedColormap


# Expt info
disk = "O:"
subject = "AF4"  # Animal
recordSession = ["004", "005"] # Unit
recordPlane = ["000", "001"]
hueExptId = ["004_008", "005_007"]  # Stimulus test
oriExptId = ["004_005", "005_010"]

pValue = 0.05

oriaucThres = 0.7
diraucThres = 0.7
hueaucThres = 0.9

mainpath = joinpath(disk,subject, "2P_analysis")
dataFolder = joinpath(mainpath, join(["U", oriExptId[1][1:3]]), join([oriExptId[1], "_", recordPlane[1]]), "DataExport")
dataExportFolder1 = joinpath(mainpath, "Summary", "DataExport")
roibkgFile=matchfile(Regex("[A-Za-z0-9]*[A-Za-z0-9]*_[A-Za-z0-9]*_roibkg.jld2"),dir=dataExportFolder1,adddir=true)[1]
roi=load(roibkgFile,"roi")[1]
bkg=load(roibkgFile,"bkg")[1]

polygons=spl.geometry.polygon(roi[1])




 allROIs[i] = Polygon([tuple(l) for l in list(segment[i]-1)])


ori_lut = ListedColormap(sio.loadmat(os.path.join(colormappath, "ori_lut_alpha0.mat"), squeeze_me=True, struct_as_record=False)["lut"])
# cmap_patch = plt.cm.get_cmap("hsv")
pycall(pylt.cm.get_cmap,"hot")
cmap_patch_cpi = pylt.cm.get_cmap("jet")
cmap_patch_osi = plt.cm.get_cmap("jet")
cmap_patch_hue = plt.cm.get_cmap("hsv_r")
cmap_patch_hue1 = ListedColormap(hdf5storage.loadmat(os.path.join(colormappath,"dkllut.mat"), squeeze_me=True, struct_as_record=False)["lut"])
# if colorSpace  == "DKL":
color_hex_keys = plcolor.ListedColormap(["g", "r", "c", "m", "y", "k", "w"],name="from_list")
color_hex_keys = plcolor.ListedColormap(["#FF8080", "#FF80BF", "#FF80FF", "#BF80FF", "#8080FF", "#80BFFF", "#80FFFF",
         "#80FFBF", "#80FF80", "#BFFF80", "#FFFF80", "#FFBF80", "#808080"])   # DKL hues
# elif colorSpace == "CIE":
    color_hex_keys = plcolor.ListedColormap(["#af1600", "#8a4600", "#5a5d01", "#2a6600", "#006a00", "#006931", "#006464",
             "#0058b6", "#002DFF", "#6a2ade", "#97209b", "#aa1c50", "#808080"])   # HSL hues
# gray = "#B4B4B4"    # "#808080""#808080"

hueList = [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
fig, ax = pylt.subplots()

# if white_BK
ax.imshow(whitBK, cmap="gray", alpha=BK_alpha, vmin=0, vmax=1)
file_name_hue = os.path.join("%s_whtBK%s" % (file_name_hue, str(white_BK)))
# end
# else:
    ax.imshow(align, cmap="gray")
for vi in np.arange(numCell)
    if ((hueax_p[vi] >auc_thres) | (huedir_p[vi]>auc_thres)) & visResp[vi]:
       # ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=cmap_patch_hue(opt_hue[vi]/360)))
       ax.add_patch(PolygonPatch(allROIs[vi], alpha=1, color=color_hex_keys(hueList.index(max_hue[vi]))))
   end
end

aa=((5, 6), (10, 14), (15, 10), (10, 5), (5, 5))

ax.add_patch(desca.PolygonPatch(aa, alpha=1))
# ax.set_title( file_name_hue, fontsize=font_size)
ax.set_rasterized(true)
plt.axis("off")

pylt.savefig("visualResponsive_hist.png", dpi=300, bbox_inches="tight", pad_inches=0, format="png")



fig = figure()
art3d = PyObject(PyPlot.art3D)

# I use semicolons instead of commas, becaues I think of these as column
# vectors.  It doesn't really matter, but I guess I'm a purist.
xc = [0;0;1;1]
yc = [0;1;1;0]
zc = [1;1;2;2]
verts = (collect(zip(xc,yc,zc)),)
p3c = PyObject(art3d.Poly3DCollection(verts))

# This is necessary in order to get to a 3D axis
ax = gca(projection="3d")
pycall(ax.add_collection3d, PyAny, p3c)

# Here I simply list them as points rather than x-y-z arrays
# Then you don't have to zip() them, and this is arguably more natural
# Just including this variant in case it's useful to anyone
p1 = [0;0;1]
p2 = [0;1;1]
p3 = [0;1;2]
p4 = [0;0;2]
verts2 = ([tuple(p1...); tuple(p2...); tuple(p3...); tuple(p4...)],)

p3c2 = PyObject(art3d.Poly3DCollection(verts2, alpha=0.5))
face_color = [1, 0, 0]
pycall(p3c2.set_facecolor, PyAny, face_color)
pycall(ax.add_collection3d, PyAny, p3c2)
ax.view_init(50,30)
zlim(0,2)
