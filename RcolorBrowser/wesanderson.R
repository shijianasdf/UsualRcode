install.packages("wesanderson")
library(wesanderson)
library(help="wesanderson")
#BottleRocket1, BottleRocket2, Rushmore1, Royal1, Royal2, Zissou1, Darjeeling1, Darjeeling2,
#Chevalier1 , FantasticFox1 , Moonrise1, Moonrise2, Moonrise3, Cavalcanti1, GrandBudapest1, 
#GrandBudapest2, IsleofDogs1, IsleofDogs2
wes_palette("Cavalcanti1", 3)
wes_palette("Royal1")

# simple barplot
barplot(c(2,5,7), col=wes_palette(n=3, name="GrandBudapest2"))
