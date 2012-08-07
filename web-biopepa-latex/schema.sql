drop table if exists entries;
create table entries (
  id integer primary key autoincrement,
  title string not null,
  text string not null,
  latex string
);

drop table if exists users;
create table users (
  id integer primary key autoincrement,
  name string not null,
  password string not null
);
