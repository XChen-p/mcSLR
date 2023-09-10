function zin = zcorr(zin, mask)

    if nargin == 1
        mask = true(size(zin));
    end
    mix =   ggmfit(reshape(zin(mask),1,[])/std(zin(mask)), 3, 'ggm', [], 10)
    zin =   (zin/std(zin(mask)) - mix.mus(1))/mix.sig(1);
