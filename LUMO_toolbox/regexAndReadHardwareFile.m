function toml_data = regexAndReadHardwareFile(filename)
%This function is meant to hardware toml files so that they are
%readable by matlab-toml.
%
%Input should be the same as toml.read
%
%This functionality shouldn't be required, but matlab-toml has issues with
%reading certain valid toml files.

raw_text = fileread(filename);

%Removes indentation within toml files.
fixed_raw_text = regexprep(raw_text,'[\n\r]+[\t ]+','\n');

toml_data = toml.decode(fixed_raw_text);
end