# fitModel2Data

Fits models to data. 

# Usage

If your model conforms to the form `[a,b,c...] = function(p,x)` where p is a structure array containing parameters, and your data contains fields called `response` and `stimulus`, you can fit your model to data using:

```matlab
p = fitModel2Data(@model,data);

```

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
