function nse = NSE(obs, sim)
sse = sum((obs - sim).^2);
ssu = sum((obs - mean(obs)).^2);
nse = 1 - sse/ssu;
end