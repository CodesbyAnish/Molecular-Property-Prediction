#IMPORTING LIBRARIES AND DATA FILES
library('dplyr')
library('xgboost')
library('tidyverse')
library('Metrics')
library('ggplot2')
library('gridExtra')
library('grid')
library('corrplot')
train <- read.csv("C:/Users/anish/Downloads/Predictive Analytics/Final project/champs-scalar-coupling/train.csv") 
test <- read.csv("C:/Users/anish/Downloads/Predictive Analytics/Final project/champs-scalar-coupling/test.csv" )
dipole <- read.csv("C:/Users/anish/Downloads/Predictive Analytics/Final project/champs-scalar-coupling/dipole_moments.csv")
structure <- read.csv("C:/Users/anish/Downloads/Predictive Analytics/Final project/champs-scalar-coupling/structures.csv")
potential <- read.csv("C:/Users/anish/Downloads/Predictive Analytics/Final project/champs-scalar-coupling/potential_energy.csv")
muliken <- read.csv("C:/Users/anish/Downloads/Predictive Analytics/Final project/champs-scalar-coupling/mulliken_charges.csv")

#TRAINING DATA FRAME PROPERTIES
summary(train)
str(train)

# DATA PREPARATION
sc <- train$scalar_coupling_constant
train$scalar_coupling_constant <- NULL
full <- rbind(train,test) %>% 
  left_join(structure, by = c("molecule_name","atom_index_0" = "atom_index")) %>% 
  left_join(structure, by = c("molecule_name","atom_index_1" = "atom_index")) %>%
  mutate(
    x_dist = x.x - x.y, 
    y_dist = y.x - y.y,  
    z_dist = z.x - z.y, 
    dist = sqrt(x_dist^2 + y_dist^2 + z_dist^2)) 
names(full) <- c("id", "molecule_name","atom_index_0", "atom_index_1", "joint_type","atom_0","x_0", "y_0", "z_0", "atom_1",
                 "x_1","y_1","z_1", "x_dist","y_dist", "z_dist","distance")
train_f <- full[1:nrow(train),c("molecule_name","atom_index_0","atom_index_1","joint_type","distance")]
train_f$scalar_coupling_constant <- as.vector(sc)

#EDA 1 : FILES 
summary(train_f)
summary(test)
summary(structure)

#EDA 2 SCALAR COUPLING CONSTANT
ggplot(train_f, aes(x = scalar_coupling_constant)) +
  geom_histogram(color = "black", fill = "orange") + ggtitle("Frequency Histogram of Scalar Coupling Constant")


#EDA 3 ORIENTATION DISTRIBUTION
s1 <- nrow(subset(train_f,train_f$scalar_coupling_constant<0))
s2 <- nrow(subset(train_f,train_f$scalar_coupling_constant>0))
df <- data.frame(Orientation=c("Opposite spin","Same spin"),
                 Molecule_count=c(s2,s1))
ggplot(data=df, aes(x=Orientation, y=Molecule_count,fill=Orientation)) +
  geom_bar(stat="identity") + ggtitle("Atomic Orientation Distribution") 


#EDA 4 : CORRELATION PLOT
train1 <- full[1:nrow(train),c("molecule_name","atom_index_0","atom_index_1","joint_type",
                               "distance")]
train1$scalar_coupling_constant <- as.numeric(sc)
train2 <- train1[, c("atom_index_0", "atom_index_1","distance", "scalar_coupling_constant")]
dipole$diople_moment <- sqrt(dipole$X^2 + dipole$Y^2 + dipole$Z^2) 
tr_dip <- train1 %>%
  left_join(dipole, by = c("molecule_name"))%>%
  left_join(potential, by = "molecule_name")%>%
  left_join(muliken, by = c("molecule_name", "atom_index_0" = "atom_index"))%>%
  left_join(muliken, by = c("molecule_name", "atom_index_1" = "atom_index"))
names(tr_dip) <-c ("molecule_name","atom_index_0","atom_index_1", "joint_type","distance", 
                   "coupling_const", "X", "Y", "Z", "dipole", 
                   "potential", "mulliken0", "mulliken1")
c <- tr_dip[, c("atom_index_0", "atom_index_1", "distance", "coupling_const",
                "dipole", "potential", "mulliken0", "mulliken1")]
d <- na.omit(c) 
M <- cor(d)
corrplot.mixed(M, lower.col ="black", number.cex=0.9, upper = "circle")


#EDA 5 : Scalar Coupling Constant Distribution by Atomic Type
p1 <- ggplot(filter(train_f, scalar_coupling_constant > 0,joint_type == "1JHC"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "skyblue",
            linetype = "solid") +
  labs(x="Type = 1JHC")
p2 <- ggplot(filter(train_f,scalar_coupling_constant > 0, joint_type == "2JHC"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "lightgreen",
            linetype = "solid") +
  labs(x="Type = 2JHC")
p3 <- ggplot(filter(train_f,scalar_coupling_constant > 0, joint_type == "3JHC"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "orange",
            linetype = "solid") +
  labs(x="Type = 3JHC")
p4 <- ggplot(filter(train_f, scalar_coupling_constant > 0,joint_type == "2JHH"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "red",
            linetype = "solid") +
  labs(x="Type = 2JHH")
p5 <- ggplot(filter(train_f,scalar_coupling_constant > 0, joint_type == "3JHH"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "yellow",
            linetype = "solid") +
  labs(x="Type = 3JHH")
p6 <- ggplot(filter(train_f,scalar_coupling_constant > 0, joint_type == "3JHN"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "purple",
            linetype = "solid") +
  labs(x="Type = 3JHN")
grid.arrange(p1,p2,p3,p4,p5,p6, ncol=3,
             top = textGrob("Distribution of Coupling Constant by Type for Opposite Spin",gp=gpar(fontsize=15)))
p7 <- ggplot(filter(train_f, scalar_coupling_constant < 0,joint_type == "1JHC"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "skyblue",
            linetype = "solid") +
  labs(x="Type = 1JHC")
p8 <- ggplot(filter(train_f,scalar_coupling_constant < 0, joint_type == "2JHC"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "lightgreen",
            linetype = "solid") +
  labs(x="Type = 2JHC")
p9 <- ggplot(filter(train_f,scalar_coupling_constant < 0, joint_type == "3JHC"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "orange",
            linetype = "solid") +
  labs(x="Type = 3JHC")
p10 <- ggplot(filter(train_f, scalar_coupling_constant < 0,joint_type == "2JHH"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "red",
            linetype = "solid") +
  labs(x="Type = 2JHH")
p11 <- ggplot(filter(train_f,scalar_coupling_constant < 0, joint_type == "3JHH"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "yellow",
            linetype = "solid") +
  labs(x="Type = 3JHH")
p12 <- ggplot(filter(train_f,scalar_coupling_constant < 0, joint_type == "3JHN"), aes(x = scalar_coupling_constant)) +
  geom_area(stat = "bin", 
            binwidth = 0.5, 
            colour = "black",
            fill = "purple",
            linetype = "solid") +
  labs(x="Type = 3JHN")
grid.arrange(p7,p8,p9,p10,p11,p12, ncol=3,
             top = textGrob("Distribution of Coupling Constant by Type for Parallel Spin",gp=gpar(fontsize=15)))



#Model 1: Linear Regression 
reg <- lm(scalar_coupling_constant~ distance + joint_type , data = train_f)
summary(reg)
train_f$pred <- predict(reg, newdata = train_f)
test1 <- full[(nrow(train)+1):nrow(full),]
test1$predlm <- predict(reg, newdata = test1)
prediction_lm <- test1$predlm


#Model 2: Gradient Boosting 
train1 <- xgb.DMatrix(data = as.matrix(train_f[, c("distance")]), label = train_f$scalar_coupling_constant)
test2 <- xgb.DMatrix(data = as.matrix(test1[, c("distance")]))
xgb_params <- list("objective" = "reg:linear",
                   "eval_metric" = "mae")
model <- xgb.train(params = xgb_params,
                   data = train1,
                   eta = 0.5,
                   nrounds = 100,
                   max_depth = 2,
                   subsample = 0.9,
                   colsample_bytree = 1)
train_pred_xgb <- predict(model,  newdata = train1)
prediction_xgb <- predict(model, newdata = test2)


#Error Calculations
mae <- function(error)
{
  mean(abs(error))
}
mae(train_f$pred-train_f$scalar_coupling_constant)
mae(train_pred_xgb-train_f$scalar_coupling_constant)


#Final Predictions
final <- test1$molecule_name
final <- cbind(final,prediction_lm)
final <- cbind(final,prediction_xgb)
head(final)



