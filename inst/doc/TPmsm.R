## ----message=FALSE------------------------------------------------------------
library("TPmsm");
setThreadsTP(2);
seed <- c(2718, 3141, 5436, 6282, 8154, 9423);
setPackageSeedTP(seed);
sim_data_exp <- dgpTP(n = 1000, corr = 0, dist = "exponential",
  dist.par = c(1, 1), model.cens = "uniform", cens.par = 3,
  state2.prob = 0.5);

## -----------------------------------------------------------------------------
transAJ(object = sim_data_exp, s = 0.5108, t = 0.9163, conf = TRUE,
  conf.level = 0.95, n.boot = 1000);
transPAJ(object = sim_data_exp, s = 0.5108, t = 0.9163, conf = TRUE,
  conf.level = 0.95, n.boot = 1000);

## -----------------------------------------------------------------------------
setPackageSeedTP(seed);
sim_data_exp2 <- dgpTP(n = 1000, corr = 1, dist = "exponential",
  dist.par = c(1, 1), model.cens = "uniform", cens.par = 3,
  state2.prob = 0.5);

## -----------------------------------------------------------------------------
transKMW(object = sim_data_exp2, s = 0.5108, t = 0.9163, conf = TRUE,
  conf.level = 0.95, n.boot = 1000);
transKMPW(object = sim_data_exp2, s = 0.5108, t = 0.9163, conf = TRUE,
  conf.level = 0.95, n.boot = 1000);

## -----------------------------------------------------------------------------
data("colonTP", package = "TPmsm");
head( head(colonTP[ , c(1:4, 7)]) );

## -----------------------------------------------------------------------------
colon_obj <- with( colonTP, survTP(time1, event1, Stime, event, age) );
colon_obj_TP <- transKMW(object = colon_obj, s = 365, t = 1096,
  conf = TRUE, conf.level = 0.95);
colon_obj_TP;
colon_obj2_TP <- transKMPW(object = colon_obj, s = 365, t = 1096,
  conf = TRUE, conf.level = 0.95);
colon_obj2_TP;

## ----label=Figure2, echo=FALSE, out.height='4in', out.width='4in'-------------
colon_obj_TP <- transKMW(object = colon_obj, s = 365, conf = TRUE,
  conf.level = 0.95);
plot(colon_obj_TP, col = seq_len(5), lty = 1, ylab = "p_hj(365,t)");

## ----eval=FALSE---------------------------------------------------------------
#  colon_obj_TP <- transKMW(object = colon_obj, s = 365, conf = TRUE,
#    conf.level = 0.95);
#  plot(colon_obj_TP, col = seq_len(5), lty = 1, ylab = "p_hj(365,t)");

## ----label=Code3, eval=FALSE--------------------------------------------------
#  plot(colon_obj_TP,  tr.choice = "1 2", conf.int = TRUE, ylim = c(0, 0.2),
#    legend = FALSE, ylab = "p12(365,t)");

## ----label=Figure3, echo=FALSE, out.height='4in', out.width='4in'-------------
plot(colon_obj_TP,  tr.choice = "1 2", conf.int = TRUE, ylim = c(0, 0.2),
  legend = FALSE, ylab = "p12(365,t)");

## ----label=Code4, echo=FALSE--------------------------------------------------
CTP_obj <- transIPCW(colon_obj, s = 365, t = 1096, x = c(40, 68),
  conf = TRUE, n.boot = 1000, method.boot = "percentile");

## ----label=Figure4, echo=FALSE, out.height='4in', out.width='4in'-------------
plot(CTP_obj, plot.type = "c", tr.choice = "1 1", conf.int = TRUE,
  xlab = "Age", legend = FALSE, ylab = "p11(365,1096|age)");

## ----label=Figure5, echo=FALSE, out.height='4in', out.width='4in'-------------
plot(CTP_obj, plot.type = "c", tr.choice = "1 2", conf.int = TRUE,
  xlab = "Age", legend = FALSE, ylab = "p12(365,1096|age)");

## ----fig.keep='none'----------------------------------------------------------
CTP_obj <- transIPCW(colon_obj, s = 365, t = 1096, x = c(40, 68),
  conf = TRUE, n.boot = 1000, method.boot = "percentile");
CTP_obj;
plot(CTP_obj, plot.type = "c", tr.choice = "1 1", conf.int = TRUE,
  xlab = "Age", legend = FALSE, ylab = "p11(365,1096|age)");
plot(CTP_obj, plot.type = "c", tr.choice = "1 2", conf.int = TRUE,
  xlab = "Age", legend = FALSE, ylab = "p12(365,1096|age)");

## ----label=Code6, eval=FALSE--------------------------------------------------
#  plot(CTP_obj, plot.type = "c", col = seq_len(5), lty = 1, xlab = "Age",
#    ylab = "p_hj(365,1096|age)");

## ----label=Figure6, echo=FALSE, out.height='4in', out.width='4in'-------------
plot(CTP_obj, plot.type = "c", col = seq_len(5), lty = 1, xlab = "Age",
  ylab = "p_hj(365,1096|age)");

## -----------------------------------------------------------------------------
data("bladderTP", package = "TPmsm");
head(bladderTP);

## ----label=Code7, echo=FALSE--------------------------------------------------
bladderTP_obj <- with( bladderTP, survTP(time1, event1, Stime, event) );

## ----label=Code8, echo=FALSE--------------------------------------------------
LS2_obj <- transLS(object = bladderTP_obj, s = 3, t = 60, h = c(0.0001, 1),
  nh = 100, ncv = 100, conf = TRUE);

## ----label=Figure7, echo=FALSE, out.height='4in', out.width='4in'-------------
plot(LS2_obj, col = seq_len(5), lty = 1, ylab = "p_hj(3,t)");

## ----label=Figure8, echo=FALSE, out.height='4in', out.width='4in'-------------
plot(LS2_obj, tr.choice = "1 2", conf.int = TRUE, ylab = "p12(3,t)",
  ylim = c(0, 0.325), legend = FALSE);

## -----------------------------------------------------------------------------
bladderTP_obj <- with( bladderTP, survTP(time1, event1, Stime, event) );
LS_obj <- transLS(object = bladderTP_obj, s = 3, t = 8, h = c(0.0001, 1),
  nh = 100, ncv = 100, conf = TRUE);
LS_obj;

## ----eval=FALSE---------------------------------------------------------------
#  LS2_obj <- transLS(object = bladderTP_obj, s = 3, t = 60, h = c(0.0001, 1),
#    nh = 100, ncv = 100, conf = TRUE);
#  plot(LS2_obj, col = seq_len(5), lty = 1, ylab = "p_hj(3,t)");
#  plot(LS2_obj, tr.choice = "1 2", conf.int = TRUE, ylab = "p12(3,t)",
#    ylim = c(0, 0.325), legend = FALSE);

