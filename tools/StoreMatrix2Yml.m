function is_ok = StoreMatrix2Yml(filenameOUT, data, rowLabels, colLabels,numdigit)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% is_ok = StoreMatrix2Yml(data, rowLabels, columnLabels, numdigit )   %%%
%%%                                                                     %%%
%%%     IN: filenameOUT -->     STRING                                  %%%
%%%                             output file name                        %%%
%%%         data        -->     MATRIX                                  %%%
%%%                             with data values                        %%%
%%%         rowLables   -->     CELL ARRAY of STRINGS                   %%%
%%%                             with names of rows                      %%%    
%%%         colLables   -->     CELL ARRAY of STRINGS                   %%%
%%%                             with names of columns                   %%% 
%%%         numdigit   -->      INT                                     %%%
%%%                             number of digits                        %%% 
%%%                                                                     %%%    
%%%                                                                     %%%    
%%%     Author:     Marco Caimmi                                        %%%
%%%                 STIIMA Nationl Research Council of Italy            %%%
%%%                 marco.caimmi@stiima.cnr.it                          %%%
%%%                 marco.caimmi@gmail.com                              %%%
%%%                                                                     %%%
%%%     Year:       2021                                                %%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumRows = length(rowLabels);
NumCols = length(colLabels);

%%%% build string with row labels
RowLabel_str = 'row_label: [';
for i = 1:NumRows
    RowLabel_str = strcat(RowLabel_str, rowLabels{i});
    if i ~= NumRows
        RowLabel_str = strcat(RowLabel_str, ", ");
    end
end
    RowLabel_str =strcat( RowLabel_str,("]\n") );


%%%% build string with column labels
ColLabel_str = "col_label: [";
for i = 1:NumCols
    ColLabel_str = strcat(ColLabel_str, colLabels{i});
    if i ~= NumCols
        ColLabel_str = strcat(ColLabel_str, ", ");
    end
end
ColLabel_str =strcat( ColLabel_str,("]\n") );

%%%%% build string with matrix data
value_str = "value: [ ";
numdigitstr = strcat("%s%.", int2str(numdigit), "f");
for i = 1:NumRows
    value_str = strcat( value_str, "[" );
    for j = 1:NumCols
        value_str = sprintf(numdigitstr, value_str, data(i,j));
        if j ~= NumCols
            value_str = sprintf("%s, ", value_str);
        else
            if i ~= NumRows
                value_str = sprintf("%s], ", value_str);
            else
                value_str = sprintf("%s] ", value_str);
            end
        end
    end
end
value_str = strcat(value_str, "]\n" );

%%%%%% write string to output file
file_id = fopen(filenameOUT, "w");

fprintf( file_id, "type: \'matrix\'\n" );
fprintf( file_id, RowLabel_str );
fprintf( file_id, ColLabel_str );
fprintf( file_id, value_str );

fclose(file_id);
is_ok = true;
end