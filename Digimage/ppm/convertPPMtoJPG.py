from PIL import Image

s = "D:/warpOut/"
d = "D:/warp3/"
# for i in range(1, 10) : 
# 	im = Image.open(s + "seam_" + str(i) + ".ppm")
# 	im.save(d + "seam_" + str(i) + ".jpg")

# for i in range(1, 10) : 
# 	im = Image.open(s + "scene_0000" + str(i) + ".jpg")
# 	im.save(d + "scene_" + str(i) + ".ppm")
# for i in range(10, 100) : 
# 	im = Image.open(s + "scene_000" + str(i) + ".jpg")
# 	im.save(d + "scene_" + str(i) + ".ppm")
# for i in range(100, 300) : 
# 	im = Image.open(s + "scene_00" + str(i) + ".jpg")
# 	im.save(d + "scene_" + str(i) + ".ppm")

for i in range(1, 300) :
	im = Image.open(s + "scene_" + str(i) + ".ppm")
	im.save(d + "scene_" + str(i) + ".jpg")
