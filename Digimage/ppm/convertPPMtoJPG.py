from PIL import Image

s = "D:/carv/"
d = "D:/carvj/"
for i in range(400) : 
	im = Image.open(s + "seam_" + str(i) + ".ppm")
	im.save(d + "seam_" + str(i) + ".jpg")