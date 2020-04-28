function newStr = Reformat(strIn)
% Author: F. W. Isen
% Copyright 2009 by Morgan & Claypool
bytes = unicode2native(strIn);

badcarat = find(bytes==136)
bytes(badcarat) = 94;

badapostrophe = find(bytes==146)
bytes(badapostrophe) = 39;

newStr = char(bytes);

