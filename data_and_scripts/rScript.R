# INPUTS:
#  series: a ts object
#  ssa.resolution: The number of eigenvalues to consider when constructing SSA groups
#  arma.params: Max parameters for the residual arma as c(p,q)
#  arma.ic: Selection criteron to pass to the arma ("aic", "aicc", "bic")
#  n_ahead: number of forecasts ahead to make
#  forecast.level: CI of the forecasts to display
#  forecast.start: starting point of the forecasts; should be "series.end+1", refer to ts docs
#  arma.verbose: outputs info about arma residuals if TRUE
#  series.name: Name of the series for plotting
ssarma <- function(series, ssa.resolution, arma.params, arma.ic, n_ahead, forecast.level, forecast.start, arma.verbose, series.name){
    library("Rssa")
    library("forecast")
    ser_ssa <- ssa(series, kind="1d-ssa")
    ser_gro <- grouping.auto.wcor(ser_ssa, groups=1:ssa.resolution)
    ser_clusts <- length(ser_gro)
    ser_recon <- reconstruct(ser_ssa, groups=ser_gro)

    ser_trd <- rep(0, length(series))
    for(i in 1:length(ser_gro)){
      if(length(ser_gro[[i]] == 1)){
        ser_trd <- ser_trd + ser_recon[[i]]
      }
    }
    ser_det <- rep(0, length(series))
    for(i in 1:ser_clusts){
        ser_det <- ser_det + ser_recon[[i]]
    }
    res <- (series - ser_det)

    res_arma <- auto.arima(res, d=0, D=0, max.p=arma.params[1], max.q=arma.params[2], stationary=TRUE, seasonal=FALSE, ic=arma.ic, allowmean=FALSE)

    if(arma.verbose){
        print(res_arma)
        checkresiduals(res_arma, plot=FALSE)
    }

    for_ssa = rforecast(ser_ssa, groups = ser_gro, len = n_ahead, only.new = TRUE)
    for_det <- rep(0, n_ahead)
    for(i in 1:ser_clusts){
        for_det = for_det + for_ssa[[i]]
    }

    for_arma = forecast(res_arma, h=n_ahead, level=95)

    for_syn_pt = ts(for_det + for_arma$mean, start=forecast.start, frequency=frequency(series))
    for_syn_lw = ts(for_det + for_arma$lower, start=forecast.start, frequency=frequency(series))
    for_syn_up = ts(for_det + for_arma$upper, start=forecast.start, frequency=frequency(series))

    for_t = time(for_syn_lw)
    plot.ts(cbind(series, for_syn_pt, for_syn_lw, for_syn_up), ylab="series", plot.type="single", col=c("black", "blue", "darkblue", "darkblue"), main=sprintf("%s Series %d-month Forecasting, %d CI", series.name, n_ahead, forecast.level))
    polygon(c(for_t, rev(for_t)), c(for_syn_lw, rev(for_syn_up)), col=rgb(0,0,1,0.25), border=NA)

    list(ssa=ser_ssa, groups=ser_gro, recon=ser_recon, trend=ser_trd, arma=res_arma, forecasts=list(mean=for_syn_pt, lower=for_syn_lw, upper=for_syn_up, t=for_t), synth=ts(c(series, for_syn_pt), start=start(series), frequency=frequency(series)))
}

comparemodels <- function(data1, model1, data2, model2, y.label, series.name){
  plot.ts(cbind(data1, model1$forecasts$mean, data2, model2$forecasts$mean), main=sprintf("%s Comparison", series.name), ylab=y.label, col=c("black", "red", "black", "blue"), plot.type="single")
  polygon(c(model1$forecasts$t, rev(model1$forecasts$t)), c(model1$forecasts$lower, rev(model1$forecasts$upper)), col=rgb(1,0,0,0.25), border=NA)
  polygon(c(model2$forecasts$t, rev(model2$forecasts$t)), c(model2$forecasts$lower, rev(model2$forecasts$upper)), col=rgb(0,0,1,0.25), border=NA)
}

sigcon.comparemodels <- function(data1, model1, data2, model2, sigcon, t.inter, y.label, series.name){
  plot.ts(cbind(data1, model1$forecasts$mean, data2, sigcon$mean), main=sprintf("%s Comparison", series.name), ylab=y.label, col=c("black", "red", "black", "purple"), plot.type="single")
  abline(v=t.inter, lty=3)
  polygon(c(model1$forecasts$t, rev(model1$forecasts$t)), c(model1$forecasts$lower, rev(model1$forecasts$upper)), col=rgb(1,0,0,0.25), border=NA)
  polygon(c(time(sigcon$lower), rev(time(sigcon$lower))), c(sigcon$lower, rev(sigcon$upper)), col=rgb(1,0,1,0.25), border=NA)
  legend("bottomleft", legend=c("data", "pre-COVID model", "post-COVID model"),col=c("black", "red", "purple"), lwd=2, cex=1.2)
}

# Pred vectors must have same length: w+1:p
# Both sigma2 must be constant
#
sigmoid.converge <- function(true.p, true.s, approach.p, approach.s, start.t, pred.f){
  pred.means <- sigmoid.converge.mean(true.p, approach.p, start.t, pred.f)
  pred.devs <- 1.96 * sigmoid.converge.sd(true.s, approach.s, start.t, pred.f, length(true.p))
  pred.lower <- pred.means - pred.devs
  pred.upper <- pred.means + pred.devs
  list(mean=pred.means, lower=pred.lower, upper=pred.upper)
}

sigmoid.converge.mean <- function(true.preds, approach.preds, start.time, pred.freq){
  len <- length(true.preds)
  q <- rep(0, len)
  for(i in 1:length(q)){
    q[i] <- (sig.hat(i, len) * true.preds[i]) + ((1-sig.hat(i, len)) * approach.preds[i])
  }
  ts(q, start=start.time, frequency=pred.freq)
}

sigmoid.converge.sd <- function(true.sigma2, approach.sigma2, start.time, pred.freq, len){
  v <- rep(0, len)
  
  for(i in 1:len){
    v[i] <- sqrt(((sig.hat(i, len)^2) * true.sigma2) + (((1-sig.hat(i, len))^2) * approach.sigma2))
  }
  ts(v, start=start.time, frequency=pred.freq)
}

sig.hat <- function(x, l){
  xp <- (1/3)*(x-(l/2))
  1/(1+exp(1)^(-xp))
}


#####################
### PROJECT USAGE ###
#####################


data = read.csv("092023AdjustedSumsT.csv")
bus_ts = ts(data["Bus"], start=2002, frequency=12)
rail_ts = ts(data["Rail"], start=2002, frequency=12)
ferry_ts = ts(data["Ferry"], start=2002, frequency=12)
other_ts = ts(data["Other"], start=2002, frequency=12)

rail_ts = rail_ts + other_ts

b_bus_ts = ts(bus_ts[1:218], start=2002, frequency=12)
b_rail_ts = ts(rail_ts[1:210], start=2002, frequency=12)
b_fer_ts = ts(ferry_ts[1:218], start=2002, frequency=12)

b_bus_model <- ssarma(b_bus_ts, 23, c(5,5), "aicc", 96, 95, c(2020, 3), FALSE, "Before-times Bus")
b_rail_model <- ssarma(b_rail_ts, 24, c(5,5), "aicc", 120, 95, c(2019, 7), FALSE, "Before-times Rail")
b_fer_model <- ssarma(b_fer_ts, 22, c(5,5), "aicc", 96, 95, c(2020, 3), FALSE, "Before-times Ferry")
checkresiduals(b_bus_model$arma)
checkresiduals(b_rail_model$arma)
checkresiduals(b_fer_model$arma)

a_bus_ts = ts(bus_ts[234:length(bus_ts)], start=c(2021,6), frequency=12)
a_rail_ts = ts(rail_ts[234:length(rail_ts)], start=c(2021,6), frequency=12)
a_fer_ts = ts(ferry_ts[234:length(ferry_ts)], start=c(2021,6), frequency=12)

a_bus_model <- ssarma(a_bus_ts, 5, c(5,5), "aicc", 34, 95, c(2023, 10), FALSE, "After-times Bus")
a_rail_model <- ssarma(a_rail_ts, 5, c(5,5), "aicc", 52, 95, c(2023, 10), FALSE, "After-times Rail")
a_fer_model <- ssarma(a_fer_ts, 5, c(5,5), "aicc", 45, 95, c(2023, 10), FALSE, "After-times Ferry")
checkresiduals(a_bus_model$arma)
checkresiduals(a_rail_model$arma)
#checkresiduals(a_fer_model$arma)

bus_sigcon <- sigmoid.converge(b_bus_model$forecasts$mean[44:77], b_bus_model$arma$sigma2, a_bus_model$forecasts$mean[1:34], a_bus_model$arma$sigma2, c(2023, 10), 12)
rail_sigcon <- sigmoid.converge(b_rail_model$forecasts$mean[53:103], b_rail_model$arma$sigma2, a_rail_model$forecasts$mean[1:52], a_rail_model$arma$sigma2, c(2023, 10), 12)
fer_sigcon <- sigmoid.converge(b_fer_model$forecasts$mean[44:89], b_fer_model$arma$sigma2, a_fer_model$forecasts$mean[1:45], a_fer_model$arma$sigma2, c(2023, 10), 12)

sigcon.comparemodels(b_bus_ts, b_bus_model, a_bus_ts, a_bus_model, bus_sigcon, (2026+(6/12)), "ridership", "Bus")
sigcon.comparemodels(b_rail_ts, b_rail_model, a_rail_ts, a_rail_model, rail_sigcon, (2027+(12/12)), "ridership", "Rail")
sigcon.comparemodels(b_fer_ts, b_fer_model, a_fer_ts, a_fer_model, fer_sigcon, (2027+(5/12)), "ridership", "Ferry")


plot.ts(cbind(ts(b_bus_model$forecasts$mean[1:43], start=c(2020,3), frequency=12), bus_sigcon$mean, ts(b_bus_model$forecasts$mean[78:96], start=c(2026,7), frequency=12)), col=c("red", "purple", "red"), plot.type="single")





