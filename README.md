# fitModel2Data

Fits models to data. This is a wrapped for MATLAB's constrained optimization algorithms. 

# Usage

If your model conforms to the form `[R,a,b,c...z] = function(S,p)` where `p` is a structure array containing parameters, and `R` is a vector and `S` is a matrix, and your data contains vector fields called `response` and `stimulus`, you can fit your model to data using:

```matlab
p = fitModel2Data(@model,data);
```

and `fitModel2Data` will fit your model to your data using the [patternsearch](https://www.mathworks.com/help/gads/patternsearch.html) algorithm. 

You can also use [fmincon](https://www.mathworks.com/help/optim/ug/fmincon.html) to fit your model:

```matlab
p = fitModel2Data(@model,data,'engine','fmincon');
```

`fitModel2Data` is internally cached, using [cache.m](https://github.com/sg-s/srinivas.gs_mtools/blob/master/docs/cache.md) so, it remembers the best solution and can automatically start from there. 


# Installation 

`fitModel2Data` is written in MATLAB. It should work on any OS, but has only been tested on macOS. 

The best way to install `fitModel2Data` is through my package manager: 

```matlab
urlwrite('http://srinivas.gs/install.m','install.m'); 
install sg-s/fitModel2Data
install sg-s/srinivas.gs_mtools  
```

# License

[GPL v3](http://gplv3.fsf.org/)
