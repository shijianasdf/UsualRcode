install.packages("colortools")
library(colortools)
# The function wheel() can be used to generate a color wheel for a given color:
wheel("darkblue", num = 12)
#[1] "#00008B" "#46008B" "#8B008B" "#8B0046" "#8B0000" "#8B4500" "#8B8B00" "#468B00" "#008B00" "#008B45"
#[11] "#008B8B" "#00468B"

#Analogous color scheme
analogous("darkblue")

#Complementary color scheme
complementary("steelblue")

#Split Complementary Color Scheme
splitComp("steelblue")

#Tetradic Color Scheme
tetradic("steelblue")

#Square color scheme
square("steelblue")

#Sequential colors
sequential("steelblue")
