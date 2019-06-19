rm(list=ls(all=TRUE))
library(dplyr)
library(orientlib)
library(spdep)
setwd("C:/Users/Owner/Dropbox/Sijun_Research/BTW analysis")

dataSet <- read.csv("ZDEVICEMOTION_SUNDAY06092019.csv", header = T)
dataSet <- dataSet %>% mutate(gravity = sqrt(ZGRAVITYX ^ 2 + ZGRAVITYY ^2 + ZGRAVITYZ ^ 2))
# dataSet <- dataSet %>% mutate(user_acc = sqrt(ZUSERACCELX ^ 2 + ZUSERACCELY ^ 2 + ZUSERACCELZ ^ 2))

# calcualte the unit vector at the gravity direction
dataSet <- dataSet %>% mutate(unit_ZGRAVITYX = ZGRAVITYX / gravity)
dataSet <- dataSet %>% mutate(unit_ZGRAVITYY = ZGRAVITYY / gravity)
dataSet <- dataSet %>% mutate(unit_ZGRAVITYZ = ZGRAVITYZ / gravity)

# calculate the magnitude of user accelaration to the gravity 
dataSet <- dataSet %>% mutate(user_acc_to_gravity = -ZUSERACCELX * unit_ZGRAVITYX - ZUSERACCELY * unit_ZGRAVITYY - ZUSERACCELZ * unit_ZGRAVITYZ)


####################################################################################################################
# rotate matrix
dirCols <- c("ZQUATX", "ZQUATY", "ZQUATZ","ZQUATW")
rotateMatrics <- rotmatrix(quaternion(as.matrix(dataSet[, dirCols])))

dataSet$row_nums <- 1:dim(dataSet)[1]


rotateXYZ <- function(data, rotateMatrics) {
	index <- data[1]
 	x <- data[2]
	y <- data[3]
	z <- data[4]
	
	rotated_matrix <- rotateMatrics[[index]]
	inverse_rotateMatrix <- solve(rotated_matrix)
	rotated_xyz <- inverse_rotateMatrix %*% c(x, y, z)
	return(rotated_xyz)
}







getDirection <- function(angle) {
	if (angle > 180) {
		angle <- angle - 360
	}
	return(angle)
}

calculateAngle <- function(data) {
	x1 <- data[1]
	y1 <- data[2]
	x2 <- data[3]
	y2 <- data[4]

	theta <- acos((x1 * x2 + y1 * y2) / sqrt(x1 ^ 2 + y1 ^2))
	return(theta)
}

dirCols <- c("row_nums", "ZUSERACCELX", "ZUSERACCELY","ZUSERACCELZ")
outcomes <- t(apply(dataSet[, dirCols], 1, rotateXYZ, rotateMatrics = rotateMatrics))
outcomes <- data.frame(outcomes)
names(outcomes) <- c("original_x", "original_y", "original_z")
dataSet <- data.frame(dataSet, outcomes)


dirCols <- c("row_nums", "ZGRAVITYX", "ZGRAVITYY","ZGRAVITYZ")
magnifieldXYZ <- t(apply(dataSet[, dirCols], 1, rotateXYZ, rotateMatrics = rotateMatrics))
magnifieldXYZ <- data.frame(magnifieldXYZ)
names(magnifieldXYZ) <- c("x", "y", "z")
dataSet <- data.frame(dataSet, magnifieldXYZ)
sign(dataSet$ZHEADING - 180)



dataSet$time <- as.numeric(substr(dataSet$ZTIMESTAMPEST, 28,29))

dataSet$diff_z <- dataSet$original_z - dataSet$user_acc_to_gravity
summary(dataSet$diff_z)


non_duplicated_index <- which(!duplicated(dataSet$time))
plot(dataSet$original_y, type="b", xaxt = "n", ylim=c(-0.6, 0.6))
axis(1, at = non_duplicated_index, labels = dataSet$time[non_duplicated_index])
# points(sign(dataSet$magnifield_x), type="l")
abline(h = -0.3, lty = 2, col = 2)
abline(h = 0.3, lty = 2, col = 2)


plot(dataSet$original_x , type="b", xaxt = "n", ylim=c(-0.6, 0.6))
axis(1, at = non_duplicated_index, labels = dataSet$time[non_duplicated_index])
# points(sign(dataSet$magnifield_x), type="l")
abline(h = -0.3, lty = 2, col = 2)
abline(h = 0.3, lty = 2, col = 2)

dataSet$XY <- sqrt(dataSet$original_x ^ 2 + dataSet$original_y ^ 2)
plot(dataSet$XY, type="b", xaxt = "n")
axis(1, at = non_duplicated_index, labels = dataSet$time[non_duplicated_index])
abline(h = -0.3, lty = 2, col = 2)
abline(h = 0.3, lty = 2, col = 2)


non_duplicated_index <- which(!duplicated(dataSet$time))
plot(dataSet$original_y, type="b", xaxt = "n", col = 1, lty = 1, lwd = 1, ylim = c(-0.6, 0.6))
points(dataSet$original_x, type="b", col = 2, lty = 2, lwd = 2, pch = 2)
axis(1, at = non_duplicated_index, labels = dataSet$time[non_duplicated_index])
abline(h = -0.3, lty = 3, col = 3)
abline(h = 0.3, lty = 3, col = 3)



zlocations <- read.csv(file = "ZLOCATIONS_SUNDAY06092019.csv", header = T)
zlocations <- filter(zlocations, ZCOURSE != -1) # remove values equal to -1
# interploate zcourse for the dataSet
dataSet$zcourse_interpolate <- approx(zlocations$ZTIMESTAMPEPOCH, zlocations$ZCOURSE, xout = dataSet$ZTIMESTAMPEPOCH)$y
dataSet <- dataSet[!is.na(dataSet$zcourse_interpolate), ]

vehicle_dimensions <- Rotation(dataSet[, c("original_x", "original_y")], dataSet$zcourse_interpolate * pi / 180)
vehicle_dimensions <- data.frame(vehicle_dimensions)
names(vehicle_dimensions) <- c("vehicle_x", "vehicle_y")
dataSet <- data.frame(dataSet, vehicle_dimensions)


pdf("accelerations.pdf", height = 9, width = 12)
    plot(dataSet$vehicle_x, type="b", xaxt = "n", ylim=c(-1, 1))
    axis(1, at = non_duplicated_index, labels = dataSet$time[non_duplicated_index])
    abline(h = -0.3, lty = 2, col = 2)
    abline(h = 0.3, lty = 2, col = 2)

    plot(dataSet$vehicle_y, type="b", xaxt = "n", ylim=c(-1, 1))
    axis(1, at = non_duplicated_index, labels = dataSet$time[non_duplicated_index])
    abline(h = -0.3, lty = 2, col = 2)
    abline(h = 0.3, lty = 2, col = 2)
dev.off()

k1 <- dataSet$original_x ^ 2 + dataSet$original_y ^ 2 
k2 <- dataSet$vehicle_x ^ 2 + dataSet$vehicle_y ^ 2
sum(abs(k1 - k2))/length(k1)











dataSet$zcourse_interpolate_to_xDegree <- (360 - dataSet$zcourse_interpolate + 90) %% 360

dataSet$zcourse_x <- cos(dataSet$zcourse_interpolate_to_xDegree)
dataSet$zcourse_y <- sin(dataSet$zcourse_interpolate_to_xDegree)
dataSet$angle <- apply(dataSet[, c("original_x", "original_y", "zcourse_x", "zcourse_y")], 1, calculateAngle) * 180 / pi

hist(dataSet$angle)

plot(dataSet$angle)

n <- 200
plot(dataSet$row_nums[1:n], dataSet$calculated_new[1:n],  type = "l", col = 1, lty = 1, lwd = 1)
points(dataSet$row_nums[1:n], dataSet$zheading_new[1:n], type = "l", col = 2, lty = 1, lwd = 1)
points(dataSet$row_nums[1:n], zlocations$ZCOURSE_new[1:n], type = "l", col = 3, lty = 1, lwd = 1)







x <- 1:10
y <- rnorm(10)
par(mfrow = c(2,1))
plot(x, y, main = "approx(.) and approxfun(.)")
points(approx(x, y), col = 2, pch = "*")
points(approx(x, y, method = "constant"), col = 4, pch = "*")



coords <- cbind(1, 0)

coords.90 <- Rotation(coords, pi / 6)
coords.90





