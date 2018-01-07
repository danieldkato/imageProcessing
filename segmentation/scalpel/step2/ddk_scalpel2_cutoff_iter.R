library("scalpel")
inpt = "/mnt/nas2/homes/dan/MultiSens/data/test_movies/5036-2_test_movie/segmentation/step1out.Rdata"
min = 0.2
increment = 0.1

load(inpt)

for (i in seq(min, 1.0, increment)){
	step2out = scalpelStep2(step1Output=step1out, cutoff=i)
}
