function toml_data = regexAndReadRecordingDataFile(filename)
%This function is meant to recordingData toml files so that they are
%readable by matlab-toml.
%
%Input should be the same as toml.read
%
%This functionality shouldn't be required, but matlab-toml has issues with
%reading certain valid toml files.

raw_text = fileread(filename);

%Fixes arrays with new lines and indentation before numbers or opening
%brackets.
raw_text = regexprep(raw_text,'[\n\r]+[\t ]+([\d\[])','$1');

%Fixes arrays with newlines before closing bracket.
raw_text = regexprep(raw_text,'[\n\r]+\]',']');

%Removes spaces between delimiters and numeric elements of an array.
raw_text = regexprep(raw_text,'([^=]) +(\d+)','$1$2');

%Removes spaces between numeric elements of an array and delimters.
raw_text = regexprep(raw_text,'(\d+) +([^=])','$1$2');

%Adds spaces between, delimiters and the numeric element that comes after.
raw_text = regexprep(raw_text,',(.)',', $1');

toml_data = toml.decode(raw_text);
end