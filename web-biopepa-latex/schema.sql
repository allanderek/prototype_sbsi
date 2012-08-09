drop table if exists entries;
create table entries (
  ident integer primary key autoincrement,
  owner_id integer,
  visibility string not null,
  title string not null,
  modelsource string not null,
  latex string,
  modelsbml string 
);

drop table if exists users;
create table users (
  ident integer primary key autoincrement,
  name string not null,
  password string not null
);
