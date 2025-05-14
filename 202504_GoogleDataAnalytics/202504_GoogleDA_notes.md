**Notes for Google Data Analytics Professional Certificate**

// (Very condensed notes taken here... mostly around spreadsheet and SQL basics)



Course link: <https://www.coursera.org/professional-certificates/google-data-analytics>


---

The overall data analysis process:

> Ask -> prepare -> process -> analyze -> share -> act 

Each step corresponds to the focus of each course~~ let's begin!

---

## Course 2: Ask Questions to Make Data-Driven Decisions

Ask more questions!!!!

### Spreadsheet Basics

After we've asked the questions and defined what we need to do with the data, we'll turn to spreadsheets to help build evidence and visualize to support the findings.

#### Spreadsheet tasks:

- Organize the data

	- Pivot table

		- Sort and filter

- Calculate the data

	- Formula: a set of instructions that performs a specific calculation

		- Start the formula with an equal sign, type the expression, and Enter. Can also output Boolean values.

		- Cell reference: a cell or a range of cells in a worksheet that can be used in a formula, containing the letter of the column and the number of the row for that cell.

			- Changes made to the referenced cells will cause the formula's output to update automatically.

			- When a formula is copied to a new cell (*fill handle*: drag the tiny box in the lower right-hand corner of the cell; automatically adds the formula/function to the cells filled), its output will be updated automatically since the cells referenced in the formula will be changed accordingly.

	- Functions

		- A *preset* command. 

		- `=func()` + Enter

		- Also updated after copied to another cell

	- Errors

		- =IFERROR(formula, "text shown if error")

		- #ERROR! (only in Google Sheets) parsing error; a formula can't be interpreted as input

		- #N/A (null; an empty field)

		- #NAME? - a formula or function name isn't recognized

		- #NUM! - a formula or function calculation can't be performed as specified

		- #VALUE! - a *general* error indicating a problem with a formula or referenced cells

		- #REF! - a formula is referencing a cell that is no longer valid or has been deleted (can be avoided by referencing a range of cells (also via ':') instead of direct reference of specific cells)



---

Before you **communicate** (eg. by email, meeting, or presentation), think about

 - Who ur audience is

 - What they already know

 - What they need to know

 - How you can communicate that effectively to them

Best projects start when the communication is clear about what's expected. Knowing exactly where to start and what you need to do allows for much more efficiency.

To deal with conflicts, you can:

- Reframe the problem. Try asking, 'how can I help you reach ur goal', instead of focusing on what went wrong or who to blame. 

- Discussion is key to conflict resolution. Start a conversation or ask things like are there other important things I should be considering? (gives a chance for ur team members or stakeholders to fully lay out concerns)

- Understand the context. Focus on the end goal.


## Course 3: Prepare Data for Exploration

### Data types and structures

#### Formats of data

- **First/second/third-party** data: collected by myself/collected by another group and then sold/sold by a provider that didn't collect the data.

- **Discrete/continuous** data: quantitative data that cannot (eg. cents)/can (eg. time taken) have almost any numeric value.

- **Nominal/ordinal** data: qualitative data that is categorized *without/with a set of order or scale*.

- **Internal/external** (primary/secondary) data: generated within/outside of an organization.

- **Structured data**: organized in a certain format like *rows and columns* (eg. spreadsheets, relational databases); **unstructured** data: not organized in any easily identifiable manner (most data, eg. audio/video files, emails, photos).

#### Data tables

- Rows/columns -> '**records'/'fields**' (sometimes 'field' also refers to a single cell; a single piece of information)

- **Wide** data: every data subject has a single row with multiple columns to hold the values of various attributes of the subject. **Long** data: each row is one *time point*/observation per subject, so each subject will have data in multiple rows.

### Databases

Where data is stored.

- **Metadata**: data about data (not the data itself, but where it comes from, when/how it was created, etc. for consistency and uniformity)

	- Descriptive/structural/administrative metadata: identifies that piece of data/indicates its structure/indicates the technical source of a digital asset.

	- Stored in metadata repositories.

- **Relational database**: a database that contains a series of *related tables that can be connected via their relationships* (one or more of the *same columns* must exist inside the related tables).

- **Schema**: a way of describing how something is organized.

- **Primary key**: an identifier that references a column in which *each value for the row is unique*, cannot be null. Identifies a record in a relational database table. Only one primary key is allowed in a table. **Foreign key**: a field within a table that is a primary key in another (related) table.

	- Distinguish the **unique key** from the primary key! The former *allows for null values, and multiple unique keys are allowed in a table*, eg. both phone number and email; they're not a *unique identifier* of the row in a *relational database*. 

#### Import, sort and filter data

- CSV: comma-separated values. Use plain text.
  
	- Separations can be auto-detected by the spreadsheet app.

- Anytime you're **sorting** (arranging in meaningful order) data, it's a good idea to **freeze the header row** first (select the first row, 'freeze -> one row' from the **'View' menu**). The header row now stays visible when we scroll down the sheet.

	- Select the column -> drop-down arrow -> sort alphabetically. Or select all -> 'sort range from Col_ to Cell_' from the 'Data' menu, and click 'Data has header row'.

- **Filtering**: 'Create a filter' from the 'Data' menu, filter button will appear in the corner of each column header. Rows that don't meet the conditions are temporarily **hidden**; to make them visible again, just click the filter button to turn off.

- **Import data dynamically**: importing from the 'File' menu will only import static files and won't be updated automatically. In Google Sheets the `IMPORT` *functions* like `=IMPORTRANGE("URL", "sheet_name!cell_range")`, `IMPORTHTML` (i.e. scraping to extract data from public web pages) and `IMPORTDATA` allow us to import data on the web dynamically. 

	- `=QUERY(sheet and range, "Select *")` enables ** pseudo SQL statements** in Google Sheets to import and extract data from another spreadsheet.

- In contrast to spreadsheets, SQL does not include a function for importing data. Instead, a method we can use is the `INSERT INTO` command together with a `SELECT` statement:
```SQL
INSERT INTO destinationTable
SELECT column_name(s)
FROM sourceTable
WHERE condition
```

#### SQL and BigQuery

**SQL**: Structured Query Language. Data analysts use query languages to communicate with large amounts of data (even trillions of rows).

SQL is a language used to interact with database programs like Oracle MySQL or Microsoft SQL Server. Usually used when the data is super large or when there are multiple files within a database, since unlike spreadsheets, it can pull information from different sources in the database.

Capitalization usually doesn't matter in SQL, so `sElECt` will also work. Some SQL dialects are case-sensitive eg. BigQuery, whereas MySQL and some others aren't. Comments are added by `--` for single-line and `/* ... */` for multi-line.

**Queries**: 

1. `SELECT`: choose the *columns* (`*`: all columns). `snake_case` names for columns. 

2. `FROM`: choose the *tables* where the columns you want are located. `CamelCase` names for tables. 

3. `WHERE`: filter the *rows*  where the conditions are met. (in the example below, 'source' is a field name (a header)). **The `WHERE` clause is where we can see conditional operators** (eg. =/>/<).

```SQL
WHERE source = 'Phone'
```

Use `AND` + condition in the next line below `WHERE` to add multiple filtering conditions.

---

Just came back from the future on 2025-04-19: SQL's written order is designed for human readability, but the database engine executes the query in a different order, usually like this: **FROM** (and **JOIN**) to identify data source -> **WHERE** to filter rows from the source -> **GROUP BY** to prepare for aggregation -> **HAVING** to filter groups by rows -> **SELECT** to compute expressions and select columns -> **DISTINCT** to remove duplicates -> **ORDER BY** -> **LIMIT**

`~~/ᐠ -˕-マ ᶻ 𝗓 𐰁 ⋆。𖦹°‧ ₊˚⊹˚ ༘⋆˚`


## Course 4: Process Data from Dirty to Clean

Data integrity: the accuracy, completeness, consistency, and trustworthiness of data throughout its lifecycle.

- For internal data, usually kept by data engineers and data warehousing specialists.

- Ways to address insufficient data (limited source/keeps updating/outdated/geographically-limited): identify trends with data available, wait for more if time allows, talk with stakeholders and adjust the objective, or look for a new dataset.
  
- Usually, you need a *statistical power* of at least 0.8 or 80% (indicates 80% chance) to consider ur results statistically significant. 

### Data cleaning

Dirty data: data that is incomplete, incorrect, or irrelevant to the problem you're trying to solve. Involves text errors, inconsistent labels or data types, nulls, duplicates, etc.

- Removing unwanted data and duplicates (back up the dataset first by clicking the dataset name in the bottom bar and 'Duplicate').

- Removing extra spaces and blanks.

- Fixing misspellings/capitalizations/punctuation and other typos.

- Clearing inconsistent formats.

Data merging: combining two or more datasets into a single dataset. They should be cleaned to the same standard.

- Compatibility: how well two or more datasets are able to work together.

- Data mapping: the process of matching *fields* (and defining the desired format of the fields) from one data source to another so that they would be compatible when merging.


### Data cleaning with spreadsheet tools

- **Conditional formatting**: changes how cells appear (eg. color) when their values in the specified range (for spreadsheets, highlight/click the first column and then press *Shift + click the last column in the range to select all columns in between*. Cmd + click to select multiple discrete columns) meet specific conditions. 'Format' -> 'Conditional formatting'.

- **'Data'** menu -> 'Remove duplicates': automatically searches for and eliminates duplicate rows.

- Select the columns -> **'Format'** menu -> useful tools for standardizing formats.

- Split: use detected (or specified) *delimiters* to split a string in entries of one column to substrings in multiple columns. Select *one* column -> 'Data' menu -> 'Split text to columns'. *This can also be used to standardize a mostly-numerical column with some mismatched string values eg. `"20"` to get clean numerical entries.*

#### Some spreadsheet functions

- =COUNTIF(range of cells, "value"): returns the number of cells that match the specified value. eg. `=COUNTIF( I2:I72, ">500")`

	- Don't forget to add quotation marks around the condition!

	- `=COUNTIFS` extends the fucntion to have multiple conditions.

- =LEN(cell): counts the length of the string in that cell.

- =LEFT(cell, substring length): returns the substring with the specified length that begins from the left of the string in the cell. =RIGHT() follows the same syntax.

	- =MID(cell, reference starting point, substring length) returns the segment within the string.

- =CONCATENATE(range of cells): concat the substrings. Can be useful when transforming the data into a consistent format. 


- =TRIM(cell): get rid of the *leading, trailing, and repeated* (but not those between words) spaces in the text.

#### Viewing data

- **Pivot table**: a data summarization tool used in data processing. Select the cells (each column requires a header) -> **'Insert'** (rather than 'Data' as in the course) menu -> 'Pivot table'. 

- **=VLOOKUP()**: vertical lookup; a function that searches for a certain value in a range, matches the row, and *returns the value in the specified column*. `=VLOOKUP(value_to_search_for, 'datasheet_name'!range_of_cells, column_index, FALSE)`. 

	- The *exclamation mark indicates that we are referencing a cell in a different spreadsheet from the one we're currently working in*. 

	- The 'column_index' argument is an integer that specifies the column in the range (counted from the leftmost column of the range) *containing* the value to *return*. 

	- The 'FALSE' indicates that we're looking for an exact match. When 'TRUE', returns an approximate match (equal or greater than).

- **Plotting**. Very useful when trying to identify any *skewed data or outliers*. Select the column(s) -> 'Insert' menu -> 'Chart'.

### Cleaning data in SQL

> When you run the queries, you get the results almost instantly. But it's fascinating to see if you think conceptually how much analysis the computer is doing for you based on that little bit of command code you wrote, and it's just so powerful.

> SQL is like a pseudo-programming language, but even more fun to master.

#### Some more SQL queries

- Finding out how many rows (records) match our search criteria: combine `COUNT` and `WHERE` (equivalent to the =COUNTIF() function in spreadsheets)

- Inserting new entries into columns (`val1` corresponds to `col1`)
```SQL
INSERT INTO `dataset_name`
	(col1, col2, ...)
VALUES
	(val1, val2, ...)
```	

-  Updating the value of a specified cell (here we use unique_key to identify the row/record). WHERE executes and filters the row according to the 'unique_key = ...' condition **before** SET operates on the column.
```SQL
UPDATE `dataset_name`
SET col_name = 'value'
WHERE unique_key = ...
```
- Creating a new table for a database (running a SQL query stores the data we extract in our local memory; to save the table, we'll need to download it as a spreadsheet (.csv) or save it as a new table in the database. The latter applies to situations eg. when you want to pull a trend on a regular basis, since the table can *refresh every time with the query you've written*. That way, you can directly download the results whenever you need them for a report (reminder: a report is static while a dashboard monitors live data)).
```SQL
CREATE TABLE IF NOT EXISTS
```

- Or sometimes when you have created too many of these in the database (but do be cautious)
```SQL
DROP TABLE IF EXISTS
```

#### Queries for cleaning data

- Removing duplicates: including `DISTINCT` in the `SELECT` statement
```SQL
SELECT  -- extracts the column 'col_name' from all rows
	DISTINCT col_name  -- removes duplicate values from the extracted set
FROM
	`dataset_name`
```

- `LENGTH/LEN(column)`: double-check that our string variables in the column are consistent with the known string length (say, for 'column' it's 2 in the example here). 
```SQL
SELECT
	LENGTH(column) AS new_col
FROM
	`dataset_name`
```
`AS`: annotate the returned column (in this case, the column with each cell value being the length of the corresponding string) as 'new_col'.

- **Filtering data** using output string lengths:
```SQL
SELECT
	column
FROM
	`dataset_name`
WHERE
	LENGTH(column) > 2
```
Now in the results, we can get the entries that *violated* length = 2. But we still want to select those that contain a specific substring (even though the lengths are mismatched). 

- `SUBSTR(column, letter_to_start, substring_length) = 'contained_substring'` to get around this inconsistency problem.

	- Let's get familiar with writing the SQL query!! >w< **Start by writing the `SELECT-FROM-WHERE` basic structure -> then use the SUBSTR() function in `WHERE` and add `=` to filter the rows. (// *look look* we always have to have a *condition* in the WHERE clause since how are we gonna filter otherwise)
```SQL
SELECT 
	column
FROM
	`dataset_name`
WHERE
	SUBSTR(column,1,2) = '__'
```
- `TRIM(column) = 'trimmed_value'`: eliminate *extra* spaces, query written similarly as `SUBSTR()`.
 ```SQL
SELECT 
	column
FROM
	`dataset_name`
WHERE
	TRIM(column) = 'text_with_no_extra_spaces'
```
- `CAST(column_name AS data_type)`: convert anything from **one data type to another**. 

	- If we are trying to sort a column into **descending order**: `ORDER BY column_name DESC` after `SELECT-FROM` (if you want an ascending order, just remove the `DESC` in the query). If the data type isn't the one we expect, weird stuff might happen ('sorted alphabetically'): eg. for strings 89 shows up before 799, since *strings start sorting with the first letter* rather than comparing the numbers.

	- So we use CAST() to convert it **`as`** another column and then sort the new column (or just `SELECT` and `ORDER BY` the same `CAST(column_name AS FLOAT64)`)

```SQL
SELECT 
	CAST(column_name AS FLOAT64) AS new_column
FROM
	`dataset_name`
ORDER BY 
	new_column DESC
```

- Another example of `CAST()`: ('date' is the data type here)
 ```SQL
SELECT 
	CAST(datetimes AS date) AS date_only,   -- happens after filtering through the WHERE clause, so this CAST only processes the filtered rows to minimize computation
	price
FROM
	`dataset_name`
WHERE
	datetimes BETWEEN '2025-04-01' AND '2025-04-20'
```

- `CONCAT(col1, col2, ...) AS new_col_name` is also used in the `SELECT` statement (since applied to column(s)).
  
	- `CONCAT_WS(separator, string1, string2,...)`: concatenate two or more strings together **with a separator** between each string.

	- We can also concat strings with the `II` operator in between the strings.

- `COALESCE(the_col_to_check_first, the_col_to_pull_if_null) AS non_null_col`: when selecting a column and some values are **null**, returns **non-null values of another designated column in the same row**. Again, used in the `SELECT` statement.
```SQL
SELECT 
	COALESCE(product, product_code) AS product_info
FROM
	`dataset_name`
```

### Verify and report on cleaning results

Verification: confirming that a data-cleaning effort was well-executed and the resulting data is accurate and reliable.

- Go back to the original dirty dataset and compare it to the current one to identify the remaining problems. Eg. `FIND`.

- Take a big-picture view of ur project.

- Checking whether there are inconsistent labels for a column: insert a **pivot table** for a spreadsheet -> 'Add' the target 'Rows' -> 'Add' the corresponding 'Values' and use `COUNTA` to count the total number of *different values that appear in the target row* (**note that the `COUNT` function only counts numerical values, so for string values it would return zeros**).

- Correcting misspellings in SQL: `CASE` statement: goes through one or more conditions and returns a value as soon as a condition is met; also used inside the `SELECT` statement. (reminder: when you write the SQL query, it's easier to start with the basic structure.) In the example below assume that we already know that certain misspellings exist in that column.
```SQL
SELECT
	customer_id,
	CASE 
		WHEN first_name = 'Tnoy' THEN 'Tony'
		WHEN first_name = 'Kiraa' THEN 'Kiara'
		ELSE first_name
		END AS cleaned_name. -- selects a new column 'cleaned_name' rather than the original 'first_name'
FROM 
	`dataset_name`
```

Documentation: tracking changes involved in the data-cleaning effort.

- Changelog: a file containing a chronologically ordered list of modifications made to a project.

- In Sheets: 'File' menu -> 'Version history', or just right-click a cell and 'Show edit history'.

- In SQL: 'Query history', click the icon to the right of each previous query to bring it up to the 'Query editor'.

	- You can also leave a comment in the statement describing the reason for a change.

- While ur team can view changelogs directly, stakeholders can't and have to rely on ur **report** to know what you did.

## Course 5: Analyze Data to Answer Questions

The goal of analysis is to **identify trends and relationships within data so you can accurately answer the question you're asking**. 

- Organize/collect data and get input from others -> format and adjust (eg. sorting and filtering) -> transform data (observing relationships between data points and making calculations).


### Organizing data in spreadsheets

- '**Sort sheet**': all of the data in a spreadsheet is sorted by the ranking of **one** specific sorted column; data *across rows is kept together*. '**Sort range**' *doesn't keep the information across rows together* (which isn't desirable at all if the cells in that column are specifically related to their rows), and nothing else in the whole spreadsheet besides the specified column is rearranged.

	- Select the column X -> 'Data' menu -> Sort sheet/range by column X 

- `=SORT(range_of_cells, column_index, TRUE/FALSE)` (TRUE for ascending order and FALSE for descending), sorts rows in the range by the specified column. The parameters are kinda similar to =VLOOKUP(), with the column_index being an integer counted from the leftmost column of the range. Eg. `=SORT(A2:O20, 2, TRUE)`.

- Customize sort order using multiple conditions/columns: Select the cells -> 'Data' menu -> Sort range -> *Advanced range sorting options* -> 'Data has header row', 'Add another sort column'.

### Formatting and adjusting data

- 'Format' menu -> tools that automatically convert values to certain formats

- `=CONVERT(cell, "current_format", "new_format")` 

	- Eg. `=CONVERT(B2, "C", "F")` to convert the temperature in cell B2 from degrees Celsius to degrees Fahrenheit.

	- When adding data to tables using a formula, *go back and paste the data in **as values** (rather than still a formula expression) so they're locked in to avoid confusion.* Right-click a new column, 'Paste special' -> 'Paste values only' to get the **static** values in the column.

- Select column -> 'Data' menu -> **'Data validation'** (useful in teamwork)

	- Adding **drop-down menus** with predetermined options:  select 'List of items' option in 'Criteria', and type in the items we want in the drop-down menu.

	- Creating **custom checkboxes**: select 'Checkbox' option in 'Criteria' -> tick 'Use custom cell values'.

	- Protecting structured data and formulas: select 'Reject input' option in 'On invalid data'.


- More advanced SQL queries: (in the query below, it pulls out multiple columns... and since we're using **aggregate functions** we must use `GROUP BY` to **group together rows**, otherwise we would get an error for mixing aggregated columns and non-aggregated columns like 'usertype')

	- `GROUP BY` statement is often used with **aggregate functions** for columns (used in `SELECT` statements, `COUNT()`, `MAX()`, `MIN()`, `SUM()`, `AVG()`) to group **by one or more columns**. Returns a table thingy.

	- The query below tries to answer what the most popular routes (start-to-end station combinations) by usertype are (since sorted by descending order), and what their average duration is. 
```SQL
SELECT
	usertype,
	CONCAT(start_station_name, " to ", end_station_name)  AS  route,
	COUNT(*)  as  num_trips,
	ROUND(AVG(CAST(tripduration  AS  int64)/60),2)  AS  duration
FROM
	`bigquery-public-data.new_york_citibike.citibike_trips`
GROUP  BY
	start_station_name, end_station_name, usertype
ORDER  BY
	num_trips  DESC
LIMIT  10
```
`ROUND(AVG(CAST(tripduration  AS  int64)/60),2)  AS  duration`: averaged over 60 rows and rounded to two decimal places.

#### String functions in spreadsheets (similar to those in SQL)

- `=LEN(cell)` returns the length of the string.

- `=FIND()` to **locate** specific characters in a string (case-sensitive); returns an integer as the position of the found character (the n-*th*)

- `=LEFT/RIGHT(cell, substring_length)` to pull out the substring from either end.

- `=VALUE(cell)` converts a **text string** that represents a numerical value (date/time/number) into an actual number.

### Data Aggregation

The process of gathering data from multiple sources in order to combine it into a single summarized collection. 

- **`=VLOOKUP`**

	- **Starts at the top of the specified range and searches downward vertically in each cell until the value in the first argument is found.**

		- So if it's approximate matching (set 'TRUE'), the returned value from the column to search for is actually the same one or any searched value **greater** than the designated lookup value in the first argument,,,,, 

	- Common mistakes while applying: data type inconsistencies ('Format' menu or functions like =VALUE()); extra spaces (=TRIM()); duplicates ('Data' menu -> 'Remove duplicates' tool).

	- Often used when populating data in one spreadsheet from another (remember the aforementioned **exclamation mark** used to designate the range from *a different* spreadsheet): `=VLOOKUP(value_to_search_for, 'datasheet_name'!range_of_cells, column_index, FALSE)`

	- Limitations: only returns the first value found in the range; only returns values to the right of the column it looks up (so 'Insert 1 left' and paste the column you want to search there); fragile to column changes due to index as input rather than direct column referencing (this is a convention that dates back to early versions of spreadsheets); fragile to changes in the referenced spreadsheet (you either lock the spreadsheet by 'Data' menu -> 'Protected sheets & ranges' (might annoy team-mates) or use `MATCH`, a function used to **locate the position of a specific lookup value** to help with spreadsheet version control).

	- Absolute/relative reference. Eg.`another_dataset!$A:$B`. Determines whether the row number or column letter should stay fixed (absolute) or adjusted (relative) when the formula is **dragged or filled**.

		- Using absolute referencing for the table_array in `VLOOKUP` ensures that the **range won't shift unexpectedly when you drag the formula**. 

			- Eg. the original formula in cell C2 is `=VLOOKUP(A2, another_dataset!A$2:B$100, FALSE)`, when dragged to C3, $ ensures that the range won't shift to A3:B101, which skips row 2 and potentially misses our data (since `VLOOKUP` searches from the top of the column). 

- Using **`JOIN`** in SQL to aggregate data in databases

	- Combines *rows* from two or more tables based on *a related column*, which can be seen as the SQL version of `VLOOKUP` (both **vertically match a value and return by row**).

	- Uses the aforementioned concepts of primary keys (a unique, non-null column) and foreign keys (columns that are primary keys in other/related tables). 

	- `INNER JOIN` (default for `JOIN`, only returns records where the tables are overlapping, i.e. **have the same values for a key/column in both tables**), `LEFT/RIGHT JOIN` (returns **all** records from the left/right table and only the matching records from the right/left table; the first table appeared is 'left'), `FULL OUTER JOIN`. If there are records in one table without a match, **it will create a record with null values for the other table**.
```SQL
SELECT
	*
FROM
	tableA
LEFT JOIN
	tableB ON tableA.keyA = tableB.keyB
```

- `=COUNT` and `COUNT` in spreadsheets and SQL:

	- Former: counts the total number of *numerical values* within a specific range.

	- Latter: returns the number of **rows** in a specific range; `COUNT(DISTINCT col_name)`: a query that only returns the distinct values in the range.

- Subqueries: queries nested in larger queries. The statement containing the subquery can also be called outer query/outer select. **The inner query executes first**. Eg. *useful for separating aggregate functions and non-aggregate functions in `SELECT` statements for expressions on columns*. 

- In the query, you can use aliases `AS` in both `SELECT` and `FROM` statements (so also for tables!). But you can also just omit the `AS`, see below (yahh although the code won't work,,, just note the 's' and 't' for short table aliases)
```SQL
SELECT
  station_id,
  name,
  number_of_rides number_of_rides_starting_at_station
FROM
  (
    SELECT
      start_station_id, 
      COUNT(*) number_of_rides
    FROM
      bigquery-public-data.new_york_citibike.citibike_trips
    GROUP BY
      start_station_id
  ) t
INNER JOIN 
  bigquery-public-data.new_york_citibike.citibike_stations s
  ON s.station_id = t.start_station_id
ORDER BY
  number_of_rides DESC
```
 We get an error: 

 > No matching signature for operator = for argument types: STRING, INT64 Signature: T1 = T1 Unable to find common supertype for templated argument <T1> Input types for <T1>: {INT64, STRING} at [16:74]

Because we're trying to join two columns of different data types, and it can't match/compare a string to a number! So we first use `CAST` to troubleshoot:

```SQL
INNER JOIN 
  bigquery-public-data.new_york_citibike.citibike_stations 
  ON citibike_stations.station_id = CAST(station_num_trips.start_station_id AS STRING)
```
And then no data displayed!!!!!! So I checked the two tables and found that these two columns/keys don't match at all: while the start_station_id is a 3-digit integer, the station_id in the other table is stuff like 66dd4f31-0aca-11e7-82f6-3863bb44ef7c. So extremely weird. Seems that the databases have changed.

So, similarly, alternative queries doing the same matching won't work as well:
```SQL
SELECT station_id, name
FROM bigquery-public-data.new_york_citibike.citibike_stations
WHERE 
  station_id IN
  (
    SELECT 
      start_station_id
    FROM
      bigquery-public-data.new_york_citibike.citibike_trips
    WHERE
      usertype = 'Subscriber'
  )
```

- `HAVING`:  In one of the previous course notes (where i wrote i came back from the future...) we've known that in the actual execution, `GROUP BY` executes **after** `WHERE` has filtered out the rows in the *underlying table*/raw data (operates on individual rows and can use any column from the tables in the `FROM` clause), so `HAVING` is like doing the same thing as `WHERE` (filtering rows, in this case, **groups**) but it's **processed after `GROUP BY`**, which means it can be used with aggregate functions (eg. `COUNT()`, `MAX()`, `MIN()`, `SUM()`, `AVG()`, see below).  
```SQL
SELECT department, COUNT(*) employee_count, AVG(salary) avg_salary
FROM employees_dataset
WHERE hire_date > '2020-01-01'  -- Filters individual rows first, executes before GROUP BY
GROUP BY department
HAVING COUNT(*) > 5  -- Filters the **groups** with more than 5 records/employees, processed after GROUP BY but before SELECT
ORDER BY avg_salary DESC;
```
- `CASE`: returns records after matching our conditions by allowing us to include `WHEN`, `THEN`, `AND`, `END` statements in the query. If 'ELSE' is omitted, non-matching statuses become NULL.

	- Processed left-to-right, stopping at first true condition (short-circuiting). 

	- Can appear in `WHERE`, `HAVING`, `SELECT`, and `ORDER BY` clauses.
```SQL
SELECT 
    employee_id,
    name,
    department,
    salary
FROM 
    employees
WHERE 
    CASE 
        WHEN department = 'Engineering' THEN salary >= 90000
        WHEN department = 'Sales' THEN salary >= 60000
        ELSE salary >= 45000  -- default threshold for all other departments
    END
ORDER BY 
    -- Prioritize sorting by departments, then sort high-to-low salary ('Engineering' (1) -> 'Sales' (2) -> other departments (3); within each group, sorted by salary in descending order)
    CASE 
        WHEN department = 'Engineering' THEN 1
        WHEN department = 'Sales' THEN 2
        ELSE 3
    END,
    salary DESC;
```

> ...just why can't these languages be as straightforward as pandas/Python 

(weird,, i had this complaint back in April while keeping the notes,, but rn when i came back in May to sort them i don't see the reason for it...?)

### Data calculations

#### In spreadsheets

- `=SUMIF(range_to_match_condition, 'condition', range_to_sum)` a function that adds numeric data based on one condition.  When a value in the range of the first argument meets the condition, the value in the **corresponding row** in the range of the third argument will be added to the sum. `=SUMIFS` extends the fucntion to have multiple conditions.

- `=SUMPRODUCT`: multiplies corresponding values in arrays and returns the sum of those products (kinda like vector dot products -> scalar output).

- Create a pivot table -> Add 'Rows' -> right-click a cell in the pivot table (in this example the rows are different dates) and 'Create pivot date group' -> Add 'Values', choose 'Summarize by' and 'Show as' for calculations automatically done to the groups. 

	- Creating a filter in the pivot table: add 'Filters' -> select 'Filter by condition' in 'Status'. Then set the 'Values' -> select 'COUNT'/'COUNTA' in 'Summarize by' to see how the counts differed between the rows/groups.

	- **Calculated Field**: a new field within a pivot table that carries out certain calculations (through 'Formula') based on the values of other fields, but since we've added a filter to the pivot table it's only going to calculate the range on the raw data after filtering. add 'Values' -> 'Calculated Field' -> select 'Custom' in 'Summarize by'.

#### In SQL

- Arithmetic operators: eg. `SELECT (col_A + col_B) * col_C AS new_col`

	- Mod: uses `%`. 

- Functions: similar to spreadsheets, eg. =SUM -> SUM; =AVERAGE -> AVG (these functions are aggregate functions since they **compute values on one or more input but return a single value**.)

- The notation `<>` is equivalent to `!=`. 

- `EXTRACT(datetime_part FROM col_with_timestamp)` extracts part of a timestamp value. Eg.
```SQL
SELECT
  EXTRACT(YEAR FROM starttime) AS year,
  count(*) as num_of_rides
FROM
  bigquery-public-data.new_york_citibike.citibike_trips
GROUP BY
  year
ORDER BY
  year -- default is ascending order
```
> Q: GROUP BY is based on the selected column 'year'; but shouldn't `SELECT` execute later than `GROUP BY`?

> > A: BigQuery (and some other database servers like MySQL) allows using `SELECT` aliases in `GROUP BY` as a convenience feature/*syntactic sugar*, just does extra work to map aliases to their original expressions. This syntax technically violates the standard logical execution order. 

- How `WHERE col_n != col_n` works in the data validation process: for valid data, any value in a row compared to itself should be equal so usually this won't return any value (since they would all return FALSE and get filtered out).

	- For NaN (**'Not a Number', a special floating-point value in IEEE 754 (used in SQL for numeric columns)**) among numeric values, `NaN != NaN` -> **`TRUE`** (by definition).

	- Note that `NaN` **is NOT** `NULL`, the latter means **missing data** (eg. optional input). `NULL != NULL`  -> `NULL`  (not  `TRUE`  or  `FALSE`), but  `WHERE`  only keeps rows where the condition is  `TRUE`, so  `NULL`values are excluded from results and require further detections like `WHERE col_n IS NULL` or `IS NOT NULL` (likewise, NaN also has detection methods like `IS_NAN`)

- Using **temporary tables** during calculations for convenience. The `WITH` clause is a type of temporary table* that you can query from multiple times.
```SQL
WITH trips_over_1_hr AS (
	SELECT *
	FROM bigquery-public-data.new_york_citibike.citibike_trips
	WHERE tripduration >= 60)
	
## Count how many trips are 60+ minutes long 
-- Use hashtags (one # if only for this session; two ##s (called global temporary table) allow the temporary table to be stored when you/other users start a new runtime) to denote a temporary table 

SELECT
	COUNT(*)  AS  c
FROM  trips_over_1_hr
```
- `SELECT col INTO temp_table_name FROM table_name`, another way to create temporary tables. It copies a table to another one without rlly adding the new table to the database (BigQuery doesn't recognize this command tho)

- **`CREATE TABLE` this statement DOES add the table into the database**: `CREATE TABLE table_path AS (SELECT * FROM another_table WHERE 'condition')`; note that the **table_path** must be completed to *project_id.dataset_id.table_name*. If we are using `CREATE TEMP TABLE`, a simple table_name would work and the full table path isn't needed. 

## Course 6: Share Data Through the Art of Visualization

- Dynamic vs. static visualizations.

- Having an interactive visualization can be useful for both you and the audience you share it with. But remember that the more power you give the user (more interactivity), the less control you have over the story you want the data to tell. 

- Elements: line, shape, color, space, movement. Use a headline, a subtitle and labels (while direct labeling keeps ur audience's attention on the graphic itself, legends force the audience to do more work because they are positioned further away from the chart's data). Simplify!

- Ppl in the world are visual learners. Choosing the right viz: *which one will make it easiest for the user to understand the point (5-second rule) you're trying to make?*

	- Design thinking: a process used to solve complex problems in a *user-centric* way. A design process: emphasize, define, ideate, prototype, test.

### A Tableau tour 

We're using Tableau Public in the browser. 

- Drags a lot.

- Features: data types are represented by: eg. '#' for numeric data, globe symbol for geographic data, calendar (with a clock) for date (and time) data, etc. Just click on that icon to change the data type. 

- Data storytelling: (i just realized that the annual reports on eg. NetEase Cloud Music are just extremely successful data storytelling modes-) The first step is to engage ur audience. What do they hope to get from the data insights I deliver?

	- Spotlighting: scanning through data to quickly identify the most important insights. 

	- *"Don't let the way you create something influence what it's actually saying."*

#### Tableau dashboards

**Dashboard**: a tool that organizes information from **multiple datasets into one central location** for tracking, analysis, and simple visualization. It monitors **live** data. 
 
- Vertical or horizontal layout.

- Dashboard filter: a tool for showing only the data that meets a specific criteria while hiding the rest. Filter by selecting the range and choose 'Keep only' or 'Exclude'.

### Presentations and slideshows

- Title, subtitle, presented_by, date (**it's good to specify a date eg. 'date created' or 'last modified'** to give context).

- Less than five lines and 25 words per slide. Choose ur words carefully. 

- When you include visuals on a slide, don't share too many details at once (don't make it an *eyesore chart*!!). Ask urself: *what's the **single** most important thing i want my audience to learn from my analysis?* (if there are several points, don't cram them all on one slide; instead, **create a new visual for each point**, then **add an arrow, a call-out, or another clearly-labelled element to direct** ur audience's attention toward what you want them to look at.)

	- If you **copy and paste** a visual to ur slideshow, the original file might be updated but not the visual. So **link** it, so that the visual lives within its original file, and the slideshow connects to it with the visual's URL -> **always up-to-date**. An **embedded** object is also stored in the original file, but it's a copy that **won't be updated** accordingly. 

- When you get to ur conclusion (i.e., ur big reveal/'aha' moment), ur visuals must communicate these messages with clarity AND EXCITEMENT!!! 

- Present with a framework (**a story/logical flow**). It helps keep you focused on the most important information during the presentation. 

	- Purpose statement -> tell ur story with data -> conclusion -> appendix. 

	- Start with ur understanding of the task; presenting the context will make the audience easier to connect with the data. Establish ur hypothesis early in the presentation.

	- Showcasing the metrics you used to help the audience understand the impact ur findings will have. 

	- Steps of the McCandless Method for data visualizations (moves from the general to the specific).
  
![The McCandless Method](https://img.picui.cn/free/2025/04/20/6804a93516d74.png)

- For a very big presentation, try to have an *ally* in the room who has listened to it in advance, so you can get feedback and get someone nodding and aligning to the numbers that you're about to present. 

- Don't assume that ur audience is already familiar with jargon and other stuff. Always be ready to explain further if asked. If necessary, include the explanations in the slide description. 

## Course 7: Data Analysis with R Programming

Notes kept separately in the .rmd file. 
