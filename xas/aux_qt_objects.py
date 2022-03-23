
def update_list_decorator(method):
    def wrapper(obj, *args, emit_signal=True, **kwargs):
        result = method(obj, *args, **kwargs)
        if emit_signal:
            obj.update_list_action()
        return result
    return wrapper


class PersistentList:
    list_update_signal = None

    def __init__(self, list_type=list, json_file_path='', boot_fresh=True):
        self.items = list_type()
        self.json_file_path = json_file_path
        self.init_from_settings(boot_fresh=boot_fresh)

    # @property
    # def local_file_default_path(self):
    #     return f"{ROOT_PATH}/{USER_FILEPATH}/{RE.md['year']}/{RE.md['cycle']}/{RE.md['PROPOSAL']}/"

    def init_from_settings(self, boot_fresh=False):
        if not boot_fresh:
            try:
                self.add_items_from_file(self.json_file_path)
            except FileNotFoundError:
                self.save_to_settings()

    @update_list_decorator
    def add_items_from_file(self, file):
        self.items += self.item_list_from_file(file)

    def item_list_from_file(self, file):
        with open(file, 'r') as f:
            item_list = json.loads(f.read())
        return item_list

    def save_to_settings(self):
        self.save_to_file(self.json_file_path)

    def save_to_file(self, file):
        if file:
            with open(file, 'w') as f:
                json.dump(self.items, f)

    @update_list_decorator
    def reset(self):
        self.items = []

    @update_list_decorator
    def insert(self, index, item):
        self.items.insert(index, item)

    @update_list_decorator
    def append(self, item):
        self.items.append(item)

    @update_list_decorator
    def extend(self, item_list):
        self.items.extend(item_list)

    @update_list_decorator
    def pop(self, index):
        self.items.pop(index)

    @update_list_decorator
    def pop_many(self, index_list):
        index_list.sort(reverse=True)
        for index in index_list:
            self.pop(index, emit_signal=False)

    # @update_list_decorator
    # def update_item_at_index(self, index, item):
    #     self.items[index] = index

    def __getitem__(self, index):
        return self.items[index]

    @update_list_decorator
    def __setitem__(self, index, item):
        self.items[index] = item

    def __iter__(self):
        return self.items.__iter__()

    def __repr__(self):
        return self.items.__repr__()

    def __len__(self):
        return len(self.items)

    # def item_at_index(self, index):
    #     return self.items[index]

    # def append_list_update_signal(self, signal):
    #     self.list_update_signal = signal

    def update_list_action(self):
        # if self.list_update_signal is not None:
        #     self.list_update_signal.emit()
        self.save_to_settings()

# bla = PersistentList()
# bla.append('a')
# bla[0]






# def validate_item_decorator(method):
#     def wrapper(obj, *args, emit_signal=True, **kwargs):
#         result = method(obj, *args, **kwargs)
#         if emit_signal:
#             obj.update_list_action()
#         return result
#     return wrapper

class NamedDict(dict):

    def extend(self, item_list):
        if type(item_list) != ListOfNamedDicts:
            raise Exception ('wrong type!')
        for item in item_list:
            _validate_item(item)

        if 'element_list' not in self.keys():
            self['element_list'] = []
        self['element_list'].extend(item_list)

    def append(self, item):
        _validate_item(item)
        if 'element_list' not in self.keys():
            self['element_list'] = []
        self['element_list'].append(item)

DEFAULT_KEYS = ['name']
def _validate_item(item, required_keys=DEFAULT_KEYS):
    if type(item) != dict:
        raise Exception(f'item is not a dict')
    item_keys = list(item.keys())
    if not all([req_key in item_keys for req_key in required_keys]):
        raise Exception(f'item is missing a "name" key')
    if 'element_list' in item_keys:
        for element in item['element_list']:
            validate_item(element)
    return NamedDict(item)



class ListOfNamedDicts(list):

    def __init__(self, *items):
        # if type(items) != list:
        #     items = [items]
        for item in items:
            if type(item) == list:
                for _item in item:
                    self.validate_item(_item)
                super().__init__(*items)
            else:
                self.validate_item(item)
                super().__init__(items)

    def validate_item(self, item):
        _validate_item(item)

    def insert(self, index, item):
        self.validate_item(item)
        super().insert(index, item)

    def append(self, item):
        self.validate_item(item)
        super().append(item)

    def extend(self, item_list):
        for item in item_list:
            self.validate_item(item)
        super().extend(item_list)

    def __setitem__(self, index, item):
        self.validate_item(item)
        self.items[index] = item

# bla = ListOfNamedDicts([{'name' : 'a'}])

class PersistentListWithQTreeWidget(PersistentList):

    def __init__(self, json_file_path='', boot_fresh=True):
        super().__init__(list_type=ListOfNamedDicts, json_file_path=json_file_path, boot_fresh=boot_fresh)
        self.qtwidget = None

# bla = PersistentListWithQTreeWidget()
# bla.append({'name': 'a'})
    # def append_qtwidget(self, qtwidget):
    #     self.qtwidget = qtwidget
    #
    # def _make_qtitem(self, parent, item_name, force_unchecked=False, checkable=True):
    #     qtitem = QtWidgets.QTreeWidgetItem(parent)
    #     qtitem.setText(0, item_name)
    #     qtitem.setExpanded(True)
    #     if checkable:
    #         qtitem.setFlags(qtitem.flags() | Qt.ItemIsTristate | Qt.ItemIsUserCheckable)
    #     if force_unchecked:
    #         qtitem.setCheckState(0, Qt.Unchecked)
    #     return qtitem
    #
    #
    # def insert(self, index, item):
    #     super().insert(index, item)
    #
    # def append(self, item):
    #     super().append(item)
    #
    # def extend(self, item_list):
    #     super().extend(index, item)
    #
    # def __setitem__(self, index, item):
    #     self.items[index] = item
    #
    # def update_list_action(self):
    #
    #     super().update_list_action()


